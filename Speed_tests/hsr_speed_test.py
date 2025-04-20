# Script to assess the speed performance of the hsr in constructing
# fingerprints and computing similarity scores.

from __future__ import annotations
import gzip, io, random, time, threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool, cpu_count, get_context
from typing import List, Tuple

import requests, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from tqdm import tqdm
import hsr

SDF_URL_LIST_FILE = "zinc_sdf_links.txt"
OUTPUT_CSV        = "hsr_fps_speed.csv"
TARGET_MOLECULES  = 1_000_000
MAX_THREADS       = 16
MAX_PROCESSES     = min(cpu_count(), 8)
REQUEST_TIMEOUT   = 20
RDLogger.DisableLog("rdApp.*")

stop_event = threading.Event()                    

# ── helpers ──────────────────────────────────────────────────────────
def download_mols(url: str, timeout: int = REQUEST_TIMEOUT):
    if stop_event.is_set():
        return []
    try:
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        if stop_event.is_set():
            return []
        with gzip.GzipFile(fileobj=io.BytesIO(r.content), mode="rb") as fh:
            suppl = Chem.ForwardSDMolSupplier(fh)
            mols  = [m for m in suppl if m is not None]
        return [(m, url) for m in mols]
    except Exception:
        return []

def fingerprint_mol(payload: Tuple[Chem.Mol, str]):
    mol, src = payload
    try:
        smi = Chem.MolToSmiles(mol)
        t0  = time.perf_counter()             
        fp  = hsr.fingerprint.generate_fingerprint_from_molecule(mol)
        return fp, time.perf_counter() - t0, smi, src
    except Exception:
        return None

# ── main ─────────────────────────────────────────────────────────────
def main():
    global_start = time.time()
    print(f"🎯 Target: {TARGET_MOLECULES:,} molecules  "
          f"(threads {MAX_THREADS}, processes {MAX_PROCESSES})\n")

    urls = [u.strip() for u in Path(SDF_URL_LIST_FILE).read_text().splitlines()
            if u.strip().endswith(".sdf.gz")]
    random.shuffle(urls)
    url_iter = iter(urls)

    # ── phase 1 – download ─────────────────────────────────────────
    t0_dl = time.time()
    collected: List[Tuple[Chem.Mol, str]] = []

    pool = ThreadPoolExecutor(max_workers=MAX_THREADS)
    try:
        bar = tqdm(total=TARGET_MOLECULES, desc="📥 Downloading", dynamic_ncols=True)
        futures = {pool.submit(download_mols, next(url_iter)): None
                   for _ in range(MAX_THREADS)}

        while futures and len(collected) < TARGET_MOLECULES:
            for fut in as_completed(list(futures)):
                futures.pop(fut, None)
                collected.extend(fut.result())
                bar.update(len(fut.result()))

                if len(collected) >= TARGET_MOLECULES:
                    stop_event.set()
                    for pending in futures:
                        pending.cancel()
                    futures.clear()
                    break

                try:
                    nxt = next(url_iter)
                    futures[pool.submit(download_mols, nxt)] = None
                except StopIteration:
                    pass
        bar.close()
    finally:
        pool.shutdown(wait=False, cancel_futures=True)

    collected = collected[:TARGET_MOLECULES]
    dl_wall = time.time() - t0_dl
    print(f"✅  Finished download phase "
          f"({len(collected):,} molecules from {len({s for _, s in collected}):,} files)\n")

    # ── phase 2 – fingerprint ──────────────────────────────────────
    print("🧪 Fingerprinting…")
    t0_fp_wall = time.time()
    ctx = get_context("spawn")
    with ctx.Pool(processes=MAX_PROCESSES) as mp:
        results = list(
            tqdm(mp.imap_unordered(fingerprint_mol, collected, chunksize=100),
                 total=len(collected), desc="🧪 Fingerprinting", dynamic_ncols=True))
    fp_wall = time.time() - t0_fp_wall

    fps, times, smiles, sources = zip(*[r for r in results if r])

    # ── csv ────────────────────────────────────────────────────────
    pd.DataFrame({"Source File": sources,
                  "Fingerprint Time (s)": times,
                  "SMILES": smiles}).to_csv(OUTPUT_CSV, index=False)

    # ── stats ──────────────────────────────────────────────────────
    print("\n🔍 Benchmark Statistics")
    qfp = random.choice(fps)
    t0 = time.time()
    sims = [hsr.similarity.compute_similarity_score(qfp, fp)
            for fp in fps if not np.array_equal(fp, qfp)]
    sim_time = time.time() - t0
    sort_time = time.time() - (t0 := time.time())

    print(f"✅ Molecules processed      : {len(fps):,}")
    print(f"✅ Files downloaded         : {len(set(sources)):,}")
    print(f"Average fingerprint time   : {sum(times)/len(times):.4f} s")
    print(f"Maximum fingerprint time   : {max(times):.4f} s")
    print(f"Minimum fingerprint time   : {min(times):.4f} s")
    print(f"{len(sims):,} similarity comparisons in {sim_time:.4f} s")
    print(f"Similarity sort            in {sort_time:.4f} s\n")

    print(f"⏱️  Sequential fp time         : {sum(times):.2f} s")
    print(f"⏱️  Wall‑clock fp phase        : {fp_wall:.2f} s")
    print(f"⏱️  Wall‑clock download phase  : {dl_wall:.2f} s")
    print(f"⏱️  End‑to‑end runtime         : {time.time() - global_start:.2f} s\n")

    print(f"✅ Done! Results saved to “{OUTPUT_CSV}”.\n")

if __name__ == "__main__":
    main()

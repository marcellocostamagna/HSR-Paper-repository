# Script to save all the URLs of .sdf.gz files from the ZINC database

import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

base_url = "http://files.docking.org/3D/"
output_file = "zinc_sdf_links.txt"
max_threads = 16  

def list_links(url):
    try:
        r = requests.get(url, timeout=10)
        soup = BeautifulSoup(r.text, "html.parser")
        return [a["href"] for a in soup.find_all("a") if a.get("href")]
    except Exception:
        return []

def get_leaf_folders(top_url):
    subfolders = list_links(top_url)
    return [urljoin(top_url, f) for f in subfolders if f.endswith("/")]

def get_sdf_links_from_leaf(leaf_url):
    try:
        links = list_links(leaf_url)
        return [urljoin(leaf_url, link) for link in links if link.endswith(".sdf.gz")]
    except Exception:
        return []

# Step 1: Get top-level folders
print("🔍 Scanning top-level folders in ZINC20 /3D/...\n")
top_folders = list_links(base_url)
top_folders = [urljoin(base_url, f) for f in top_folders if f.endswith("/") and f[0].isalpha()]

# Step 2: Gather all leaf folders
leaf_folders = []
print("📂 Gathering leaf folders (2 levels deep)...")
for top_url in tqdm(top_folders, desc="Top-level folders", dynamic_ncols=True):
    leaf_folders.extend(get_leaf_folders(top_url))

# Step 3: Collect all .sdf.gz links in parallel
print(f"\n🔎 Collecting .sdf.gz links from {len(leaf_folders)} folders using {max_threads} threads...\n")
sdf_links = []

with ThreadPoolExecutor(max_threads) as executor:
    futures = {executor.submit(get_sdf_links_from_leaf, url): url for url in leaf_folders}
    for future in tqdm(as_completed(futures), total=len(futures), desc="Extracting URLs", dynamic_ncols=True):
        sdf_links.extend(future.result())

# Step 4: Save to file
with open(output_file, "w") as f:
    for link in sdf_links:
        f.write(link + "\n")

print(f"\n✅ Saved {len(sdf_links)} .sdf.gz URLs to: {output_file}")

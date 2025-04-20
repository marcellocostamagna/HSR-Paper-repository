import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

base_url = "http://files.docking.org/3D/"
max_threads = 16  # Adjust depending on your system/network

def list_links(url):
    try:
        r = requests.get(url, timeout=10)
        soup = BeautifulSoup(r.text, "html.parser")
        return [a["href"] for a in soup.find_all("a") if a.get("href")]
    except Exception:
        return []

def get_leaf_folders(top_url):
    """Get all subfolders (leaves) under a given top-level folder."""
    subfolders = list_links(top_url)
    return [urljoin(top_url, f) for f in subfolders if f.endswith("/")]

def count_sdf_files_in_leaf(leaf_url):
    """Count .sdf.gz files in a single leaf folder."""
    try:
        links = list_links(leaf_url)
        return sum(1 for l in links if l.endswith(".sdf.gz"))
    except Exception:
        return 0

# Step 1: Get top-level folders (AA/, AB/, ...)
print("🔍 Scanning top-level folders in ZINC20 /3D/...\n")
top_folders = list_links(base_url)
top_folders = [urljoin(base_url, f) for f in top_folders if f.endswith("/") and f[0].isalpha()]

# Step 2: Get all leaf folders under each top-level folder
leaf_folders = []
print("📂 Gathering leaf folders (2 levels deep)...")
for top_url in tqdm(top_folders, desc="Top-level folders", dynamic_ncols=True):
    leaf_folders.extend(get_leaf_folders(top_url))

# Step 3: Count .sdf.gz files in parallel across all leaf folders
print(f"\n🔎 Counting .sdf.gz files in {len(leaf_folders)} folders using {max_threads} threads...\n")
total_sdf_files = 0

with ThreadPoolExecutor(max_threads) as executor:
    futures = {executor.submit(count_sdf_files_in_leaf, url): url for url in leaf_folders}
    for future in tqdm(as_completed(futures), total=len(futures), desc="Counting .sdf.gz", dynamic_ncols=True):
        total_sdf_files += future.result()

print(f"\n✅ Total .sdf.gz files found: {total_sdf_files}")

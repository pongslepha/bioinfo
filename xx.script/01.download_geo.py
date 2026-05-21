import os
import urllib.request
import urllib.error
import re

# 1. Define Accession list
gse_ids = ["GSE39872", "GSE134033", "GSE37114"]

# 2. Output directory
script_dir = os.path.dirname(os.path.abspath(__file__))
dest_dir = os.path.abspath(os.path.join(script_dir, "..", "00.data", "GEO"))
os.makedirs(dest_dir, exist_ok=True)
print(f"Output directory verified: {dest_dir}")

# 3. Helper progress reporter
def report_progress(block_num, block_size, total_size):
    read_so_far = block_num * block_size
    if total_size > 0:
        percent = min(100, (read_so_far * 100) // total_size)
        print(f"\rDownloading... {percent}% ({read_so_far // (1024*1024)}MB / {total_size // (1024*1024)}MB)", end="")
    else:
        print(f"\rDownloading... {read_so_far // (1024*1024)}MB", end="")

# 4. Download function with automatic fallback
def download_geo_data(gse_id, dest_dir):
    # Extract digits to construct the GEO FTP/HTTP URL structure
    match = re.match(r"GSE(\d+)", gse_id)
    if not match:
        print(f"Invalid GSE ID pattern: {gse_id}")
        return
    
    digits = match.group(1)
    if len(digits) <= 3:
        subfolder = "GSEnnn"
    else:
        subfolder = f"GSE{digits[:-3]}nnn"
        
    suppl_dir_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{subfolder}/{gse_id}/suppl/"
    raw_tar_url = f"{suppl_dir_url}{gse_id}_RAW.tar"
    raw_tar_dest = os.path.join(dest_dir, f"{gse_id}_RAW.tar")
    
    print("\n" + "="*50)
    print(f"Processing Accession: {gse_id}")
    print("="*50)
    
    # Step A: Attempt primary _RAW.tar download
    if os.path.exists(raw_tar_dest) and os.path.getsize(raw_tar_dest) > 1000:
        print(f"Primary RAW tar for {gse_id} already exists. Skipping download.")
        download_success = True
    else:
        print(f"Attempting to download primary RAW tar from: {raw_tar_url}")
        download_success = False
        
        try:
            urllib.request.urlretrieve(raw_tar_url, raw_tar_dest, reporthook=report_progress)
            
            # Check size to confirm it is not an error HTML page saved as a tar file
            if os.path.exists(raw_tar_dest) and os.path.getsize(raw_tar_dest) > 1000:
                print(f"\nSuccessfully downloaded primary RAW tar: {os.path.basename(raw_tar_dest)}!")
                download_success = True
            else:
                if os.path.exists(raw_tar_dest):
                    os.remove(raw_tar_dest)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print("\nNotice: Primary _RAW.tar file not found (404).")
            else:
                print(f"\nHTTP Error during primary download: {e}")
            if os.path.exists(raw_tar_dest):
                os.remove(raw_tar_dest)
        except Exception as e:
            print(f"\nError during primary download: {e}")
            if os.path.exists(raw_tar_dest):
                os.remove(raw_tar_dest)

    # Step B: Fallback to scanning the /suppl/ directory for individual files
    if not download_success:
        print(f"Scanning directory for supplementary files: {suppl_dir_url}")
        try:
            with urllib.request.urlopen(suppl_dir_url) as response:
                html = response.read().decode('utf-8')
            
            # Find all href links
            links = re.findall(r'href="([^"?#]+)"', html)
            # Filter links: keep files, exclude parent directories, subfolders, or absolute links
            files_to_download = [
                link for link in links 
                if not link.endswith('/') and not link.startswith('http') and not link.startswith('/')
            ]
            # De-duplicate
            files_to_download = list(dict.fromkeys(files_to_download))
            
            if files_to_download:
                print(f"Found {len(files_to_download)} supplementary file(s) for {gse_id}:")
                for f in files_to_download:
                    print(f"  - {f}")
                
                for file_name in files_to_download:
                    file_url = f"{suppl_dir_url}{file_name}"
                    file_dest = os.path.join(dest_dir, file_name)
                    
                    if os.path.exists(file_dest) and os.path.getsize(file_dest) > 1000:
                        print(f"File {file_name} already exists. Skipping download.")
                        continue
                        
                    print(f"Downloading {file_name}...")
                    try:
                        urllib.request.urlretrieve(file_url, file_dest, reporthook=report_progress)
                        print(f"\nSuccessfully downloaded: {file_name}")
                    except Exception as e:
                        print(f"\nFailed to download {file_name}: {e}")
            else:
                print(f"No supplementary files found in the directory for {gse_id}.")
        except Exception as e:
            print(f"Error accessing supplementary directory: {e}")

# 5. Execute downloads
if __name__ == "__main__":
    for gse_id in gse_ids:
        download_geo_data(gse_id, dest_dir)

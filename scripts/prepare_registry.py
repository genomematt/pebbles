import requests
from pathlib import Path


def generate_registry():
    """Generates the registry including Human, Mouse, and Yeast."""
    registry = []

    # --- Human (Homo_sapiens) ---
    for i in range(15):
        acc_v = 26 + i
        aid = "GRCh38" if i == 0 else f"GRCh38.p{i}"
        registry.append(
            {"cat": "vertebrate_mammalian", "species": "Homo_sapiens", "id": aid, "acc": f"GCF_000001405.{acc_v}"})

    # Verified available GRCh37 directories
    grch37_available = [
        (13, "GRCh37"), (14, "GRCh37.p2"), (17, "GRCh37.p5"), (21, "GRCh37.p9"),
        (22, "GRCh37.p10"), (23, "GRCh37.p11"), (24, "GRCh37.p12"), (25, "GRCh37.p13"),
    ]
    for v, aid in grch37_available:
        registry.append(
            {"cat": "vertebrate_mammalian", "species": "Homo_sapiens", "id": aid, "acc": f"GCF_000001405.{v}"})

    # --- Mouse (Mus_musculus) ---
    registry.append(
        {"cat": "vertebrate_mammalian", "species": "Mus_musculus", "id": "GRCm39", "acc": "GCF_000001635.27"})
    for i in range(7):
        acc_v = 20 + i
        aid = "GRCm38" if i == 0 else f"GRCm38.p{i}"
        registry.append(
            {"cat": "vertebrate_mammalian", "species": "Mus_musculus", "id": aid, "acc": f"GCF_000001635.{acc_v}"})

    # --- Yeast (Saccharomyces_cerevisiae) ---
    registry.append({
        "cat": "fungi",
        "species": "Saccharomyces_cerevisiae",
        "id": "R64",
        "acc": "GCF_000146045.2"
    })

    return registry


def bundle_reports():
    """Downloads reports using category-aware RefSeq paths."""
    # Logic to find src/pebbles/data/assemblies/
    base_path = Path(__file__).parent.parent / "src" / "pebbles" / "data" / "assemblies"
    base_path.mkdir(parents=True, exist_ok=True)

    targets = generate_registry()
    base_ftp = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq"

    print(f"Syncing {len(targets)} assembly reports to {base_path}...")

    success_count = 0
    for item in targets:
        cat, species, acc, aid = item['cat'], item['species'], item['acc'], item['id']
        local_file = base_path / f"{acc}.report.txt"

        if local_file.exists():
            success_count += 1
            continue

        folder_name = f"{acc}_{aid}"
        url = f"{base_ftp}/{cat}/{species}/all_assembly_versions/{folder_name}/{folder_name}_assembly_report.txt"

        try:
            r = requests.get(url, timeout=15)
            # Handle GRCh37 p0 naming inconsistency specifically if still encountered
            if r.status_code == 404 and aid == "GRCh37":
                alt_url = url.replace(f"{acc}_GRCh37", f"{acc}_GRCh37.p0")
                r = requests.get(alt_url, timeout=15)

            if r.status_code == 200:
                local_file.write_text(r.text, encoding='utf-8')
                print(f" Bundled: {species} {aid} ({acc})")
                success_count += 1
            else:
                print(f" Failed: {aid} (Status {r.status_code})")
                print(f" URL tried: {url}")

        except Exception as e:
            print(f" Network error for {aid}: {e}")

    print(f"\nCompleted: {success_count}/{len(targets)} reports ready in bundle.")


if __name__ == "__main__":
    bundle_reports()
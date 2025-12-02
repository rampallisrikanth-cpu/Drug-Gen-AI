# utils.py
# Helper functions: parse CSV and VCF into simple {rsid: genotype} dictionary

import pandas as pd
import io

def parse_23andme_csv(uploaded_file):
    """
    Parse 23andMe-style raw file or simple CSV with columns:
    rsid, chromosome, position, genotype
    """
    try:
        df = pd.read_csv(uploaded_file, sep=None, engine="python")  # auto-detect separator
    except Exception:
        uploaded_file.seek(0)
        text = uploaded_file.read().decode("utf-8")
        df = pd.read_csv(io.StringIO(text))
    # Accept columns with names in lowercase too
    cols = [c.lower() for c in df.columns]
    if "rsid" in cols and "genotype" in cols:
        # normalize column names
        rs_col = df.columns[cols.index("rsid")]
        gt_col = df.columns[cols.index("genotype")]
        result = {}
        for _, row in df.iterrows():
            rs = str(row[rs_col]).strip()
            gt = str(row[gt_col]).strip().upper()
            if rs and rs.startswith("rs"):
                # normalize genotype (remove "/")
                gt = gt.replace("/", "").replace("|", "").replace(" ", "")
                result[rs] = gt
        return result
    else:
        # fallback: treat first row as gene:value table (prototype)
        # e.g., CYP2D6,CYP2C19,... on header and values in row1
        header = list(df.columns)
        # If first column contains 'rsid' text then cannot parse
        if any(str(h).lower().startswith("rs") for h in header):
            # try reading as two columns: rsid,genotype without headers
            result = {}
            for _, row in df.iterrows():
                if len(row) >= 2:
                    rs = str(row.iloc[0]).strip()
                    gt = str(row.iloc[1]).strip().upper()
                    if rs.startswith("rs"):
                        result[rs] = gt.replace("/", "").replace("|", "").replace(" ", "")
            return result
        # else gene-expression style
        row0 = df.iloc[0]
        result = {}
        for col in header:
            # Use col name as "GENE" but convert to synthetic rs-like id for prototype
            key = col.strip()
            val = str(row0[col]).strip()
            result[key] = val
        return result

def parse_vcf(uploaded_file):
    """
    Minimal VCF parser:
    Converts VCF to {rsid: genotype_string}
    Attempts to parse sample genotype if present (FORMAT and sample columns)
    If genotype fields missing, it uses REF/ALT to create a simple allele pair.
    """
    uploaded_file.seek(0)
    text = uploaded_file.read().decode("utf-8", errors="ignore").splitlines()
    records = {}
    header_cols = None
    for line in text:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header_cols = line.lstrip("#").split("\t")
            continue
        if not line.strip():
            continue
        parts = line.split("\t")
        if header_cols:
            rec = dict(zip(header_cols, parts))
            rsid = rec.get("ID", "").strip()
            ref = rec.get("REF", "").strip()
            alt = rec.get("ALT", "").strip()
            # if there is a sample column (columns beyond INFO)
            # typical columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
            genotype = None
            if "FORMAT" in rec and len(parts) > header_cols.index("FORMAT") + 1:
                fmt = rec["FORMAT"].split(":")
                sample = parts[header_cols.index("FORMAT") + 1]
                sample_vals = sample.split(":")
                # find GT index
                try:
                    gt_index = fmt.index("GT")
                    gt = sample_vals[gt_index]
                    # GT may be 0/1 or 1|0 etc. We map 0->REF allele, 1->ALT allele
                    alleles = []
                    for allele_token in gt.replace("|", "/").split("/"):
                        if allele_token == ".":
                            alleles.append(".")
                        elif allele_token == "0":
                            alleles.append(ref)
                        elif allele_token == "1":
                            alleles.append(alt.split(",")[0] if "," in alt else alt)
                        else:
                            # numeric index >1
                            try:
                                idx = int(allele_token)
                                alleles.append(alt.split(",")[idx-1])
                            except Exception:
                                alleles.append(".")
                    genotype = "".join(alleles)
                except ValueError:
                    genotype = None
            if genotype is None:
                # fallback: make genotype from REF and first ALT (heterozygous)
                if alt and ref:
                    # set heterozygous representation
                    a1 = ref
                    a2 = alt.split(",")[0]
                    genotype = (a1 + a2) if (len(a1) == 1 and len(a2) == 1) else f"{a1}/{a2}"
                else:
                    genotype = ""
            if rsid and rsid.startswith("rs"):
                records[rsid] = genotype.replace("/", "").replace("|", "")
            else:
                # ignore non-rs lines for prototype
                continue
        else:
            # No header columns; skip (bad vcf)
            continue
    return records

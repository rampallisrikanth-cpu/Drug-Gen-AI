# model.py
# Pharmacogenomic rules and scoring logic for prototype

# Mapping of rsID -> human-readable gene (for display)
RSID_TO_GENE = {
    "rs4244285": "CYP2C19",   # *2
    "rs12248560": "CYP2C19",  # *17
    "rs1065852": "CYP2D6",    # *10 / important
    "rs3892097": "CYP2D6",    # *4 (loss of function)
    "rs1799853": "CYP2C9",    # *2
    "rs1057910": "CYP2C9",    # *3
    "rs887829": "UGT1A1",     # expression assoc (proxy)
}

# Genotype -> phenotype dictionaries (very simplified, prototype)
# Use common genotype notations: "AA", "AG", "GG", "CT", etc.
GENOTYPE_TO_PHENOTYPE = {
    "CYP2C19": {
        # rs4244285 (A = *2 loss-of-function, G = *1)
        # rs12248560 (T = *17 increased function)
        # We'll interpret combined logic in scoring function instead of enumerating every combo
    },
    "CYP2D6": {
        # We'll interpret presence of loss-of-function rs3892097 as poor; rs1065852 hetero as intermediate
    },
    "CYP2C9": {
        # rs1799853 (C>T), rs1057910 (A>C)
    },
    "UGT1A1": {
        # rs887829 (T associated with lower expression depending on population)
    }
}

# For simplicity: scoring map per phenotype for drugs
PHENOTYPE_SCORE = {
    # example scores: Poor->low effectiveness for prodrugs; Poor->higher exposure risk for active drugs.
    "poor": 10,
    "intermediate": 45,
    "normal": 85,
    "rapid": 95,
    "ultrarapid": 99
}

# Drug definitions: which rsIDs to check and what kind of mapping logic to apply
DRUGS = {
    "Codeine": {
        "description": "Codeine is converted to morphine by CYP2D6. Poor metabolizers get little pain relief; ultrarapid may have toxicity.",
        "gene_rsids": ["rs1065852", "rs3892097"],  # CYP2D6
        "type": "prodrug_cytochrome2d6"
    },
    "Clopidogrel": {
        "description": "Clopidogrel requires activation by CYP2C19. Poor metabolizers have reduced antiplatelet effect.",
        "gene_rsids": ["rs4244285", "rs12248560"],  # CYP2C19
        "type": "prodrug_cytochrome2c19"
    },
    "Ibuprofen": {
        "description": "Ibuprofen is metabolized by CYP2C9. Reduced function increases exposure and side-effect risk.",
        "gene_rsids": ["rs1799853", "rs1057910"],
        "type": "metabolized_cytochrome2c9"
    },
    "Paracetamol": {
        "description": "Paracetamol is largely glucuronidated (UGT1A1) and sulfated; variants may change toxicity risk.",
        "gene_rsids": ["rs887829"],
        "type": "glucuronidation_ugt1a1"
    },
    "Omeprazole": {
        "description": "Omeprazole is metabolized by CYP2C19; poor metabolizers have higher exposure.",
        "gene_rsids": ["rs4244285", "rs12248560"],
        "type": "prodrug_cytochrome2c19"  # same logic site as clopidogrel but effect direction differs
    }
}

# Utilities used by the scoring logic
def call_cyp2d6_phenotype(found):
    """
    Simple decision logic for CYP2D6 using presence of rs3892097 (loss) and rs1065852 (reduced).
    found: dict of rs -> genotype (e.g., {'rs3892097':'GG','rs1065852':'AG'})
    Returns phenotype string: poor/intermediate/normal/ultrarapid
    """
    # crude rules for prototype:
    # if rs3892097 present as variant (e.g., allele indicates loss) -> poor
    # if rs1065852 homozygous variant -> poor, heterozygous -> intermediate
    # otherwise normal
    # Note: real-world genotyping is more complex (star alleles)
    gt_389 = found.get("rs3892097")
    gt_106 = found.get("rs1065852")
    # interpret: assume for rs3892097 (G = loss) (prototype assumption)
    if gt_389 is not None:
        if "G" in gt_389 and gt_389.count("G") == 2:
            return "poor"
        if "G" in gt_389 and gt_389.count("G") == 1:
            return "intermediate"
    if gt_106 is not None:
        # assume variant allele is A -> reduces function
        if gt_106.count("A") == 2:
            return "poor"
        if gt_106.count("A") == 1:
            return "intermediate"
    return "normal"

def call_cyp2c19_phenotype(found):
    """
    Simple rules for CYP2C19 using rs4244285 (A = *2 loss) and rs12248560 (T = *17 increased)
    """
    gt_2 = found.get("rs4244285")  # *2
    gt_17 = found.get("rs12248560")  # *17
    loss = 0
    gain = 0
    if gt_2:
        loss += gt_2.count("A")  # prototype mapping
    if gt_17:
        gain += gt_17.count("T")
    # prioritize loss
    if loss >= 2:
        return "poor"
    if loss == 1:
        return "intermediate"
    if gain >= 1:
        return "rapid" if gain == 1 else "ultrarapid"
    return "normal"

def call_cyp2c9_phenotype(found):
    """
    Prototype rules for CYP2C9 using rs1799853 (C>T) and rs1057910 (A>C)
    Loss alleles produce intermediate/poor
    """
    loss_count = 0
    gt_2 = found.get("rs1799853")
    gt_3 = found.get("rs1057910")
    if gt_2:
        # assume T is variant causing reduced function
        loss_count += gt_2.count("T")
    if gt_3:
        # assume C is variant causing reduced function
        loss_count += gt_3.count("C")
    if loss_count >= 2:
        return "poor"
    if loss_count == 1:
        return "intermediate"
    return "normal"

def call_ugt1a1_phenotype(found):
    """
    Very simplified: rs887829 T allele associated with lower UGT1A1 expression (prototype)
    """
    gt = found.get("rs887829")
    if not gt:
        return "normal"
    tcount = gt.count("T")
    if tcount >= 2:
        return "poor"
    if tcount == 1:
        return "intermediate"
    return "normal"

# High-level function to compute drug score
def compute_drug_scores(snp_dict):
    """
    snp_dict: {rsid: genotype_string}
    returns: {drug: {score:int, phenotype:..., explanation:...}}
    """
    results = {}
    # helper to fetch only rs relevant
    for drug, meta in DRUGS.items():
        rsids = meta["gene_rsids"]
        found = {}
        for r in rsids:
            if r in snp_dict:
                found[r] = snp_dict[r]
        phenotype = "normal"
        # call appropriate phenotype function
        typ = meta["type"]
        if typ == "prodrug_cytochrome2d6":
            phenotype = call_cyp2d6_phenotype(found)
            # For codeine (prodrug), poor -> low effect, ultrarapid -> high risk
            if phenotype == "poor":
                score = 10
            elif phenotype == "intermediate":
                score = 40
            elif phenotype == "normal":
                score = 85
            elif phenotype == "rapid":
                score = 95
            else:
                score = 99
        elif typ == "prodrug_cytochrome2c19":
            phenotype = call_cyp2c19_phenotype(found)
            # For clopidogrel specifically: poor -> low effect
            if drug == "Clopidogrel":
                if phenotype == "poor":
                    score = 20
                elif phenotype == "intermediate":
                    score = 50
                elif phenotype in ("rapid","ultrarapid"):
                    score = 90
                else:
                    score = 90
            else:
                # omeprazole: poor -> higher exposure (so "effectiveness" of acid suppression is higher)
                if phenotype == "poor":
                    score = 95
                elif phenotype == "intermediate":
                    score = 75
                else:
                    score = 60
        elif typ == "metabolized_cytochrome2c9":
            phenotype = call_cyp2c9_phenotype(found)
            # for ibuprofen: poor -> increased exposure (so more side-effects); we'll map to "effectiveness" roughly normal
            if phenotype == "poor":
                score = 70
            elif phenotype == "intermediate":
                score = 80
            else:
                score = 85
        elif typ == "glucuronidation_ugt1a1":
            phenotype = call_ugt1a1_phenotype(found)
            # for paracetamol: poor glucuronidation -> perhaps higher toxicity risk; keep effectiveness near normal
            if phenotype == "poor":
                score = 60
            elif phenotype == "intermediate":
                score = 75
            else:
                score = 88
        else:
            phenotype = "unknown"
            score = 50
        explanation = f"{meta['description']} Detected genotype markers: {found if found else 'none found in file.'} => phenotype: {phenotype}."
        results[drug] = {"score": score, "phenotype": phenotype, "explanation": explanation}
    return results

from flask import Flask, render_template, request
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import IUPACData

app = Flask(__name__)

def get_3_letter_protein(one_letter_seq):
    if not one_letter_seq:
        return "None"
    # Map 'M' to 'Met', 'F' to 'Phe', etc.
    three_letter = []
    for aa in one_letter_seq:
        # Get 3-letter code from Biopython's IUPAC data
        code = IUPACData.protein_letters_1to3.get(aa, aa).capitalize()
        three_letter.append(code)
    return "-".join(three_letter)

def analyze_dna(sequence):
    clean_seq = sequence.strip().upper()
    seq_obj = Seq(clean_seq)
    
    # 1. DNA Mass
    dna_mass = f"{molecular_weight(seq_obj, 'DNA'):.2f} Da"
    
    # 2. Translation & 3-Letter Code
    protein_seq = seq_obj.translate(to_stop=True)
    three_letter_prot = get_3_letter_protein(str(protein_seq))
    
    # 3. Protein Mass
    if len(protein_seq) > 0:
        analysed_prot = ProteinAnalysis(str(protein_seq))
        prot_mass = f"{analysed_prot.molecular_weight():.2f} Da"
    else:
        prot_mass = "0.00 Da"

    # 4. Expanded AMR Scout Logic
    amr_database = {
        "TGGTATGTGGAAGTTAGATTG": "mecA (Methicillin/Staph Resistance)",
        "TTCGGCATTTCGTC": "vanA (Vancomycin Resistance)",
        "GCTTTTGCACGAAA": "bla-TEM (Penicillin Resistance)",
        "ATCAGCAACTTATC": "tetA (Tetracycline Resistance)"
    }
    
    found_amr = "No known resistance markers found."
    for marker, description in amr_database.items():
        if marker in clean_seq:
            found_amr = description
            break

    return {
        "length": len(seq_obj),
        "dna_mass": dna_mass,
        "protein_mass": prot_mass,
        "translation": three_letter_prot,
        "amr_results": found_amr,
        "restriction_sites": ["EcoRI (G*AATTC)", "BamHI (G*GATCC)", "KpnI (GGTAC*C)"]
    }

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyzer', methods=['GET', 'POST'])
def analyzer():
    results = None
    if request.method == 'POST':
        seq = request.form.get('sequence', '')
        if seq:
            results = analyze_dna(seq)
    return render_template('analyzer.html', results=results)

# Ensure you have /amr and /primer routes defined as before
if __name__ == '__main__':
    app.run(debug=True)
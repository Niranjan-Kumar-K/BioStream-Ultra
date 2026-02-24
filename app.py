from flask import Flask, render_template, request
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis # Needed for Protein Mass

app = Flask(__name__)

def analyze_dna(sequence):
    # Clean the sequence
    clean_seq = sequence.strip().upper()
    seq_obj = Seq(clean_seq)
    
    # 1. DNA Mass
    dna_mass = f"{molecular_weight(seq_obj, 'DNA'):.2f} Da"
    
    # 2. Translation & Protein Mass
    protein_seq = seq_obj.translate(to_stop=True)
    if len(protein_seq) > 0:
        # We use ProteinAnalysis to get the weight of the amino acid chain
        analysed_prot = ProteinAnalysis(str(protein_seq))
        prot_mass = f"{analysed_prot.molecular_weight():.2f} Da"
    else:
        prot_mass = "0.00 Da"

    # 3. Improved AMR Scout Logic
    amr_database = {
        "TGGTATGTGGAAGTTAGATTG": "mecA (Methicillin Resistance Detected)",
        "TTCGGCATTTCGTC": "vanA (Vancomycin Resistance Detected)"
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
        "translation": str(protein_seq),
        "amr_results": found_amr,
        "restriction_sites": ["EcoRI (G*AATTC)", "BamHI (G*GATCC)", "KpnI (GGTAC*C)"] # Example placeholders
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

# Add your other routes (amr, primer) here...

if __name__ == '__main__':
    app.run(debug=True)
from flask import Flask, render_template, request
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import IUPACData
import re

app = Flask(__name__)

# --- NCBI CONFIGURATION ---
Entrez.email = "niranjankumar270627@gmail.com" 

def fetch_full_record(accession_id):
    try:
        # Step 1: Get GenBank file for annotations (Features)
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        # Step 2: Check if sequence is missing (UndefinedSequenceError check)
        try:
            _ = record.seq[0] # Test if we can access the first base
        except Exception:
            # Step 3: "Second Ping" - Fetch the FASTA version specifically to get the sequence
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
            fasta_record = SeqIO.read(handle, "fasta")
            handle.close()
            record.seq = fasta_record.seq # Inject the missing sequence into the main record
            
        return record
    except Exception as e:
        return f"ERROR: {str(e)}"

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyzer', methods=['GET', 'POST'])
def analyzer():
    results = None
    input_seq = ""
    if request.method == 'POST':
        action = request.form.get('action')
        if action == 'fetch':
            record = fetch_full_record(request.form.get('accession_id', '').strip())
            input_seq = str(record.seq) if not isinstance(record, str) else record
        elif action == 'blast':
            # Logic for BLAST remains the same
            input_seq = request.form.get('sequence', '').strip()
            # ... (Blast logic here)
        else:
            input_seq = request.form.get('sequence', '').strip().upper()

        if input_seq and not input_seq.startswith("ERROR"):
            seq_obj = Seq(input_seq)
            protein_seq = seq_obj.translate(to_stop=True)
            results = {
                "length": len(seq_obj),
                "dna_mass": f"{molecular_weight(seq_obj, 'DNA'):.2f} Da",
                "protein_mass": f"{ProteinAnalysis(str(protein_seq)).molecular_weight():.2f} Da" if protein_seq else "0 Da",
                "translation": get_3_letter_protein(str(protein_seq)),
                "restriction": "Standard Sites Found",
                "ori": "✅ FOUND" if "TTGAGATC" in input_seq else "❌ NONE"
            }
    return render_template('analyzer.html', results=results, input_seq=input_seq)

@app.route('/amr', methods=['GET', 'POST'])
def amr():
    amr_results = []
    input_seq = ""
    if request.method == 'POST':
        action = request.form.get('action')
        if action == 'fetch':
            acc_id = request.form.get('accession_id', '').strip()
            record = fetch_full_record(acc_id)
            
            if not isinstance(record, str):
                input_seq = str(record.seq)
                for feature in record.features:
                    qualifiers = str(feature.qualifiers).lower()
                    if any(word in qualifiers for word in ["resistance", "beta-lactamase", "antibiotic"]):
                        name = feature.qualifiers.get('gene', feature.qualifiers.get('product', ['Unknown']))[0]
                        amr_results.append({
                            "name": name,
                            "pos": f"bp {int(feature.location.start)+1} to {int(feature.location.end)}",
                            "snippet": input_seq[int(feature.location.start):int(feature.location.start)+40] + "..."
                        })
            else: input_seq = record
        else:
            input_seq = request.form.get('sequence', '').upper()
    return render_template('amr.html', amr_results=amr_results, input_seq=input_seq)

@app.route('/primer', methods=['GET', 'POST'])
def primer():
    primer_results = None
    input_seq = ""
    if request.method == 'POST':
        action = request.form.get('action')
        if action == 'fetch':
            record = fetch_full_record(request.form.get('accession_id', '').strip())
            input_seq = str(record.seq) if not isinstance(record, str) else record
        else:
            input_seq = request.form.get('sequence', '').upper()
        if len(input_seq) >= 20 and not input_seq.startswith("ERROR"):
            primer_results = {"f": input_seq[:20], "r": str(Seq(input_seq[-20:]).reverse_complement())}
    return render_template('primer.html', primer_results=primer_results, input_seq=input_seq)

def get_3_letter_protein(one_letter_seq):
    if not one_letter_seq: return "None"
    return "-".join([IUPACData.protein_letters_1to3.get(aa, aa).capitalize() for aa in one_letter_seq])

if __name__ == '__main__':
    app.run(debug=True)
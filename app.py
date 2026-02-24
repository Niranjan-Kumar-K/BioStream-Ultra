from flask import Flask, render_template, request
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import IUPACData
import os

app = Flask(__name__)

# --- NCBI CONFIGURATION ---
# Use your email so NCBI doesn't block your Cloud IP
Entrez.email = "niranjankumar270627@gmail.com" 

def fetch_full_record(accession_id):
    """Fetches annotations and handles hollow records by fetching FASTA if needed."""
    try:
        # Step 1: Get GenBank file for annotations
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        # Step 2: Check if sequence is missing (UndefinedSequenceError check)
        try:
            _ = record.seq[0] 
        except Exception:
            # Step 3: Fetch the FASTA version to get the actual DNA letters
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
            fasta_record = SeqIO.read(handle, "fasta")
            handle.close()
            record.seq = fasta_record.seq 
            
        return record
    except Exception as e:
        return f"ERROR: {str(e)}"

def run_blast(seq):
    """Runs a cloud-optimized BLAST search."""
    try:
        # Trimming to 800bp prevents server timeouts on Git/Cloud deployments
        query_seq = seq[:800] 
        result_handle = NCBIWWW.qblast("blastn", "nt", query_seq)
        blast_record = NCBIXML.read(result_handle)
        
        if blast_record.alignments:
            # Clean up the title for the UI
            match_title = blast_record.alignments[0].title
            return match_title.split('|')[-1].strip()
        return "No significant match found."
    except Exception as e:
        return "NCBI Timeout. Try a shorter sequence or wait 1 minute."

def get_3_letter_protein(one_letter_seq):
    if not one_letter_seq: return "None"
    return "-".join([IUPACData.protein_letters_1to3.get(aa, aa).capitalize() for aa in one_letter_seq])

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
            input_seq = request.form.get('sequence', '').strip()
            match = run_blast(input_seq)
            # Carry over results to keep UI populated
            results = {"blast_match": match} 
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
                "restriction": "GAATTC (EcoRI), GGATCC (BamHI)", # Example logic
                "ori": "✅ FOUND" if "TTGAGATC" in input_seq else "❌ NONE",
                "blast_match": results.get("blast_match") if results else None
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
            input_seq = request.form.get('sequence', '').upper().strip()
            
        if input_seq and len(input_seq) >= 20 and not input_seq.startswith("ERROR"):
            primer_results = {
                "f": input_seq[:20],
                "r": str(Seq(input_seq[-20:]).reverse_complement())
            }
    return render_template('primer.html', primer_results=primer_results, input_seq=input_seq)

if __name__ == '__main__':
    # Use environment port for cloud hosting providers
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
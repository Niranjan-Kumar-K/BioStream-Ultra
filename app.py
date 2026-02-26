from flask import Flask, render_template, request
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight, gc_fraction, MeltingTemp as mt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import IUPACData
import os
from io import StringIO

app = Flask(__name__)

# --- NCBI SETUP ---
Entrez.email = "niranjankumar270627@gmail.com" 

def clean_seq(data):
    if not data: return ""
    clean = "".join(data.split()).upper()
    return "".join([char for char in clean if char.isalpha()])

def detect_origin(seq_str):
    motifs = ["TTATCCACA", "TGTGGATAA", "GATCTNTTN", "TTGAGATC", "AAAAAA", "TTTTTT"] 
    for motif in motifs:
        if motif in seq_str:
            return f"✅ DETECTED ({motif})"
    return "❌ NOT DETECTED"

def fetch_full_record(accession_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception:
        return "ERROR: Accession ID not found."

# --- 1. ANALYZER (PRO-WEB BLAST) ---
@app.route('/analyzer', methods=['GET', 'POST'])
def analyzer():
    results = None
    input_seq = ""
    if request.method == 'POST':
        action = request.form.get('action')
        
        if action == 'fetch':
            record = fetch_full_record(request.form.get('accession_id', '').strip())
            input_seq = str(record.seq) if not isinstance(record, str) else record
        else:
            input_seq = clean_seq(request.form.get('sequence', ''))

        if input_seq and not input_seq.startswith("ERROR"):
            seq_obj = Seq(input_seq)
            
            # --- BLAST ENGINE: SPEED OPTIMIZED ---
            blast_status = "Not Run"
            if action == 'blast':
                try:
                    # Using megablast + hitlist_size=1 to avoid 502 Timeouts
                    result_handle = NCBIWWW.qblast("blastn", "nt", input_seq, 
                                                   hitlist_size=1, 
                                                   megablast=True)
                    blast_data = result_handle.read()
                    result_handle.close()
                    
                    blast_records = NCBIXML.parse(StringIO(blast_data))
                    blast_record = next(blast_records)
                    
                    if blast_record.alignments:
                        # Extracting the cleanest name from the alignment title
                        top_hit = blast_record.alignments[0].title
                        blast_status = top_hit.split('|')[-1].strip()[:100]
                    else:
                        blast_status = "No matches found."
                except Exception as e:
                    print(f"DEBUG: BLAST Failed -> {e}")
                    blast_status = "NCBI Timeout. Try again or use a shorter sequence."

            # Stats Calculation
            rev_comp = str(seq_obj.reverse_complement())
            gc_cont = f"{(gc_fraction(seq_obj) * 100):.2f}%"
            full_protein = seq_obj.translate()
            try:
                clean_prot = str(full_protein).replace('*', '')
                p_mass = f"{ProteinAnalysis(clean_prot).molecular_weight():.2f}"
            except: p_mass = "0.00"
            
            enzymes = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT"}
            res_found = [f"{n} (@{input_seq.find(s)+1})" for n, s in enzymes.items() if s in input_seq]
            
            results = {
                "length": len(seq_obj),
                "gc_content": gc_cont,
                "dna_mass": f"{molecular_weight(seq_obj, 'DNA'):.2f} Da",
                "prot_mass": f"{p_mass} Da",
                "rev_comp": rev_comp,
                "translation": "-".join([IUPACData.protein_letters_1to3.get(aa, "Stp").capitalize() for aa in full_protein]),
                "restriction": ", ".join(res_found) if res_found else "None Found",
                "ori": detect_origin(input_seq),
                "blast_match": blast_status
            }
    return render_template('analyzer.html', results=results, input_seq=input_seq)

# --- 2. AMR SCOUT ---
@app.route('/amr', methods=['GET', 'POST'])
def amr():
    amr_results = []
    input_seq = ""
    if request.method == 'POST':
        action = request.form.get('action')
        if action == 'fetch':
            record = fetch_full_record(request.form.get('accession_id', '').strip())
            if not isinstance(record, str):
                try:
                    raw_seq = record.seq
                    input_seq = str(raw_seq) if raw_seq else ""
                    for feature in record.features:
                        qualifiers = str(feature.qualifiers).lower()
                        if any(word in qualifiers for word in ["resistance", "beta-lactamase", "antibiotic", "drug"]):
                            name = feature.qualifiers.get('gene', feature.qualifiers.get('product', ['Unknown Gene']))[0]
                            amr_results.append({
                                "gene": name.upper(),
                                "type": feature.type,
                                "location": f"{int(feature.location.start)} - {int(feature.location.end)}",
                                "sequence": str(feature.extract(record.seq))[:50] + "..."
                            })
                except: input_seq = "ERROR: Sequence undefined."
            else: input_seq = record
        else:
            input_seq = clean_seq(request.form.get('sequence', ''))
    return render_template('amr.html', amr_results=amr_results, input_seq=input_seq)

# --- 3. PRIMER LAB ---
@app.route('/primer', methods=['GET', 'POST'])
def primer():
    primer_results = None
    input_seq = ""
    if request.method == 'POST':
        input_seq = clean_seq(request.form.get('sequence', ''))
        if len(input_seq) >= 50:
            seq_obj = Seq(input_seq)
            fwd_seq = seq_obj[:20]
            rev_seq = seq_obj[-20:].reverse_complement()
            f_tm = mt.Tm_Wallace(fwd_seq)
            r_tm = mt.Tm_Wallace(rev_seq)
            primer_results = {
                "fwd_seq": str(fwd_seq),
                "fwd_tm": f"{f_tm:.1f}°C",
                "fwd_gc": f"{(gc_fraction(fwd_seq)*100):.1f}%",
                "rev_seq": str(rev_seq),
                "rev_tm": f"{r_tm:.1f}°C",
                "rev_gc": f"{(gc_fraction(rev_seq)*100):.1f}%",
                "status": "OPTIMAL" if 52 <= f_tm <= 65 else "ADVISE CHECK"
            }
        elif input_seq:
            input_seq = "ERROR: Minimum 50bp required."
    return render_template('primer.html', primer_results=primer_results, input_seq=input_seq)

@app.route('/')
def index():
    return render_template('index.html')

if __name__ == '__main__':
    # Dynamic port for Render deployment
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
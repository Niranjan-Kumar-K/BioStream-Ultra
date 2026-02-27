from flask import Flask, render_template, request, redirect, url_for
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight, gc_fraction, MeltingTemp as mt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import IUPACData
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime
import pytz
import os
from io import StringIO

app = Flask(__name__)

# --- DATABASE CONFIG ---
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///biostream.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

def get_local_time():
    tz = pytz.timezone('Asia/Kolkata')
    return datetime.now(tz)

class AnalysisRecord(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    timestamp = db.Column(db.DateTime, default=get_local_time)
    tool_name = db.Column(db.String(50))
    sequence = db.Column(db.Text)
    result_summary = db.Column(db.Text)

with app.app_context():
    db.create_all()

# Ensure you use your own registered email for NCBI
Entrez.email = "niranjankumar270627@gmail.com" 

# --- UTILITIES ---
def clean_seq(data):
    if not data: return ""
    lines = data.splitlines()
    # Strip FASTA header if present
    clean_data = "".join(lines[1:]) if lines and lines[0].startswith('>') else "".join(lines)
    clean = "".join(clean_data.split()).upper()
    return "".join([char for char in clean if char.isalpha()])

def detect_origin(seq_str):
    motifs = ["TTATCCACA", "TGTGGATAA", "GATCTNTTN", "TTGAGATC", "AAAAAA", "TTTTTT"] 
    for motif in motifs:
        if motif in seq_str:
            return f"✅ DETECTED ({motif})"
    return "❌ NOT DETECTED"

def fetch_full_record(accession_id):
    try:
        # FASTA is lightweight and prevents 500 errors on fetch
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        raw_data = handle.read()
        handle.close()
        lines = raw_data.splitlines()
        sequence = "".join(lines[1:]).upper()
        return sequence
    except Exception:
        return "ERROR: Accession ID not found."

# --- ROUTES ---

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/vault')
def vault():
    all_records = AnalysisRecord.query.order_by(AnalysisRecord.timestamp.desc()).all()
    return render_template('vault.html', records=all_records)

@app.route('/clear_vault')
def clear_vault():
    try:
        AnalysisRecord.query.delete()
        db.session.commit()
    except:
        db.session.rollback()
    return redirect(url_for('vault'))

@app.route('/analyzer', methods=['GET', 'POST'])
def analyzer():
    results = None
    input_seq = ""
    too_long = False
    
    if request.method == 'POST':
        action = request.form.get('action')
        if action == 'fetch':
            input_seq = fetch_full_record(request.form.get('accession_id', '').strip())
        else:
            input_seq = clean_seq(request.form.get('sequence', ''))

        if input_seq and not input_seq.startswith("ERROR"):
            seq_obj = Seq(input_seq)
            
            # --- HIGH CAPACITY LIMIT (300k bp) ---
            if len(seq_obj) > 300000:
                too_long = True
            else:
                blast_status = "Not Run"
                if action == 'blast':
                    try:
                        search_seq = input_seq
                        is_sliced = False
                        # SMART SLICE: If sequence is > 10kb, take a 1500bp slice to avoid NCBI rejection
                        if len(seq_obj) > 10000:
                            mid = len(seq_obj) // 2
                            search_seq = input_seq[mid:mid+1500]
                            is_sliced = True
                        
                        result_handle = NCBIWWW.qblast("blastn", "nt", search_seq, hitlist_size=1, megablast=True)
                        blast_data = result_handle.read()
                        result_handle.close()
                        
                        blast_record = next(NCBIXML.parse(StringIO(blast_data)))
                        if blast_record.alignments:
                            match_name = blast_record.alignments[0].title.split('|')[-1].strip()[:100]
                            blast_status = f"{match_name} (via slice)" if is_sliced else match_name
                        else:
                            blast_status = "No matches found."
                    except:
                        blast_status = "NCBI Timeout/Busy."

                # --- MASS CALCULATIONS (kDa) ---
                try:
                    dna_mass_val = molecular_weight(seq_obj, 'DNA') / 1000
                    dna_display = f"{dna_mass_val:.2f} kDa"
                except: dna_display = "0.00 kDa"
                
                full_protein = seq_obj.translate()
                try:
                    clean_prot = str(full_protein).replace('*', '')
                    p_mass_val = ProteinAnalysis(clean_prot).molecular_weight() / 1000
                    p_display = f"{p_mass_val:.2f} kDa"
                except: p_display = "0.00 kDa"
                
                # --- OTHER METRICS ---
                gc_cont = f"{(gc_fraction(seq_obj) * 100):.2f}%"
                rev_comp = str(seq_obj.reverse_complement())
                enzymes = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT"}
                res_found = [f"{n} (@{input_seq.find(s)+1})" for n, s in enzymes.items() if s in input_seq]
                
                # Full Protein Translation String
                translation_str = "-".join([IUPACData.protein_letters_1to3.get(aa, "Stp").capitalize() for aa in full_protein])
                
                results = {
                    "length": len(seq_obj), 
                    "gc_content": gc_cont,
                    "dna_mass": dna_display,
                    "prot_mass": p_display, 
                    "rev_comp": rev_comp,
                    "translation": translation_str,
                    "restriction": ", ".join(res_found) if res_found else "None Found",
                    "ori": detect_origin(input_seq), 
                    "blast_match": blast_status
                }
                
                db.session.add(AnalysisRecord(tool_name="Analyzer", sequence=input_seq[:500], result_summary=f"Length: {len(seq_obj)}bp"))
                db.session.commit()
                
    return render_template('analyzer.html', results=results, input_seq=input_seq, too_long=too_long)

@app.route('/amr', methods=['GET', 'POST'])
def amr():
    amr_results = []
    input_seq = ""
    if request.method == 'POST':
        action = request.form.get('action')
        if action == 'fetch':
            acc_id = request.form.get('accession_id', '').strip()
            try:
                # gbwithparts is faster for RefSeq IDs like NG_
                handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="gbwithparts", retmode="text")
                record = next(SeqIO.parse(handle, "genbank"))
                handle.close()
                input_seq = str(record.seq)
                for feature in record.features:
                    qual_str = str(feature.qualifiers).lower()
                    if any(word in qual_str for word in ["resistance", "beta-lactamase", "mec", "antibiotic"]):
                        name = feature.qualifiers.get('gene', feature.qualifiers.get('product', ['Detected Gene']))[0]
                        amr_results.append({
                            "gene": name.upper(), "type": feature.type,
                            "location": f"{int(feature.location.start)} - {int(feature.location.end)}",
                            "sequence": str(feature.extract(record.seq))[:50] + "..."
                        })
                if not amr_results:
                    for feature in record.features:
                        if feature.type == "CDS":
                            name = feature.qualifiers.get('gene', ['Unknown CDS'])[0]
                            amr_results.append({"gene": f"{name.upper()} (Potential)", "type": "CDS", "location": f"{int(feature.location.start)} - {int(feature.location.end)}", "sequence": str(feature.extract(record.seq))[:50] + "..."})
            except Exception as e:
                input_seq = f"NCBI Fetch Error: {e}"
        else:
            input_seq = clean_seq(request.form.get('sequence', ''))
        db.session.add(AnalysisRecord(tool_name="AMR Scout", sequence=input_seq[:500], result_summary=f"Found {len(amr_results)} genes"))
        db.session.commit()
    return render_template('amr.html', amr_results=amr_results, input_seq=input_seq)

@app.route('/primer', methods=['GET', 'POST'])
def primer():
    primer_results = None
    input_seq = ""
    if request.method == 'POST':
        input_seq = clean_seq(request.form.get('sequence', ''))
        if len(input_seq) >= 50:
            seq_obj = Seq(input_seq)
            fwd_seq, rev_seq = seq_obj[:20], seq_obj[-20:].reverse_complement()
            f_tm, r_tm = mt.Tm_Wallace(fwd_seq), mt.Tm_Wallace(rev_seq)
            primer_results = {
                "fwd_seq": str(fwd_seq), "fwd_tm": f"{f_tm:.1f}°C",
                "fwd_gc": f"{(gc_fraction(fwd_seq)*100):.1f}%",
                "rev_seq": str(rev_seq), "rev_tm": f"{r_tm:.1f}°C",
                "rev_gc": f"{(gc_fraction(rev_seq)*100):.1f}%",
                "status": "OPTIMAL" if 52 <= f_tm <= 65 else "ADVISE CHECK"
            }
            db.session.add(AnalysisRecord(tool_name="Primer Lab", sequence=input_seq[:500], result_summary=f"Fwd Tm: {f_tm:.1f}C"))
            db.session.commit()
    return render_template('primer.html', primer_results=primer_results, input_seq=input_seq)

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
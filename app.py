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
import re
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

Entrez.email = "niranjankumar270627@gmail.com" 

# --- UTILITIES ---
def clean_seq(data):
    if not data: return ""
    lines = data.splitlines()
    clean_data = "".join(lines[1:]) if lines and lines[0].startswith('>') else "".join(lines)
    clean = "".join(clean_data.split()).upper()
    return "".join([char for char in clean if char.isalpha()])

def detect_map_features(seq_str):
    features = []
    # Common Restriction Sites
    sites = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "NotI": "GCGGCCGC"}
    for name, site in sites.items():
        for match in re.finditer(site, seq_str):
            features.append({"name": name, "pos": match.start(), "type": "restriction"})
    
    # Origins of Replication
    ori_motifs = ["TTATCCACA", "TGTGGATAA", "TTGAGATC"]
    for motif in ori_motifs:
        for match in re.finditer(motif, seq_str):
            features.append({"name": "ORI", "pos": match.start(), "type": "origin"})
    return features

def fetch_full_record(accession_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        raw_data = handle.read()
        handle.close()
        lines = raw_data.splitlines()
        return "".join(lines[1:]).upper()
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
            if len(seq_obj) > 300000:
                too_long = True
            else:
                blast_status = "Not Run"
                if action == 'blast':
                    try:
                        search_seq = input_seq if len(seq_obj) <= 10000 else input_seq[len(seq_obj)//2:len(seq_obj)//2+1500]
                        result_handle = NCBIWWW.qblast("blastn", "nt", search_seq, hitlist_size=1, megablast=True)
                        blast_record = next(NCBIXML.parse(StringIO(result_handle.read())))
                        blast_status = blast_record.alignments[0].title.split('|')[-1].strip()[:100] if blast_record.alignments else "No matches."
                    except: blast_status = "NCBI Busy."

                gc_cont = f"{(gc_fraction(seq_obj) * 100):.2f}%"
                full_protein = seq_obj.translate()
                translation_str = "-".join([IUPACData.protein_letters_1to3.get(aa, "Stp").capitalize() for aa in full_protein])
                
                # Mapping Logic
                map_data = detect_map_features(input_seq)

                results = {
                    "length": len(seq_obj), "gc_content": gc_cont,
                    "dna_mass": f"{(molecular_weight(seq_obj, 'DNA')/1000):.2f} kDa",
                    "prot_mass": f"{(ProteinAnalysis(str(full_protein).replace('*','')).molecular_weight()/1000):.2f} kDa", 
                    "rev_comp": str(seq_obj.reverse_complement()),
                    "translation": translation_str,
                    "map_features": map_data,
                    "blast_match": blast_status
                }
                db.session.add(AnalysisRecord(tool_name="Analyzer", sequence=input_seq, result_summary=f"{len(seq_obj)}bp | GC: {gc_cont}"))
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
                handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="gbwithparts", retmode="text")
                record = next(SeqIO.parse(handle, "genbank"))
                input_seq = str(record.seq)
                for f in record.features:
                    qual = str(f.qualifiers).lower()
                    if any(w in qual for w in ["resistance", "beta-lactamase", "antibiotic"]):
                        amr_results.append({
                            "gene": f.qualifiers.get('gene', ['Unknown'])[0].upper(),
                            "type": f.type, "location": f"{int(f.location.start)} - {int(f.location.end)}",
                            "sequence": str(f.extract(record.seq))[:50] + "..."
                        })
            except: input_seq = "Fetch Error"
        else:
            input_seq = clean_seq(request.form.get('sequence', ''))
        db.session.add(AnalysisRecord(tool_name="AMR Scout", sequence=input_seq, result_summary=f"Detected {len(amr_results)} threats"))
        db.session.commit()
    return render_template('amr.html', amr_results=amr_results, input_seq=input_seq)

@app.route('/primer', methods=['GET', 'POST'])
def primer():
    primer_results = None
    input_seq = ""
    if request.method == 'POST':
        input_seq = clean_seq(request.form.get('sequence', ''))
        if len(input_seq) >= 50:
            fwd = Seq(input_seq[:20])
            rev = Seq(input_seq[-20:]).reverse_complement()
            
            f_tm, r_tm = mt.Tm_Wallace(fwd), mt.Tm_Wallace(rev)
            f_gc, r_gc = gc_fraction(fwd)*100, gc_fraction(rev)*100
            tm_diff = abs(f_tm - r_tm)
            
            f_clamp = fwd[-1] in ['G', 'C']
            r_clamp = rev[-1] in ['G', 'C']
            
            warnings = []
            if tm_diff > 5: warnings.append(f"High Tm Difference: {tm_diff:.1f}°C")
            if not f_clamp or not r_clamp: warnings.append("Missing 3' GC Clamp")
            
            status = "OPTIMAL" if not warnings else "CAUTION"
            anneal_temp = min(f_tm, r_tm) - 5
            
            protocol = {
                "mix": [
                    {"component": "Nuclease-Free Water", "amount": "32.5 µL"},
                    {"component": "10X PCR Buffer", "amount": "5.0 µL"},
                    {"component": "dNTP Mix (10mM)", "amount": "1.0 µL"},
                    {"component": "Primers (10µM each)", "amount": "5.0 µL"},
                    {"component": "Template & Taq", "amount": "1.5 µL"}
                ],
                "cycling": [
                    {"step": "Denaturation", "temp": "95°C", "time": "30s"},
                    {"step": "Annealing", "temp": f"{anneal_temp:.1f}°C", "time": "30s"},
                    {"step": "Extension", "temp": "72°C", "time": "1m/kb"}
                ]
            }

            primer_results = {
                "fwd_seq": str(fwd), "fwd_tm": f"{f_tm:.1f}°C", "fwd_gc": f"{f_gc:.1f}%", "fwd_clamp": "YES" if f_clamp else "NO",
                "rev_seq": str(rev), "rev_tm": f"{r_tm:.1f}°C", "rev_gc": f"{r_gc:.1f}%", "rev_clamp": "YES" if r_clamp else "NO",
                "tm_diff": f"{tm_diff:.1f}°C", "status": status, "warnings": warnings, "protocol": protocol
            }
            db.session.add(AnalysisRecord(tool_name="Primer Lab", sequence=input_seq, result_summary=f"F:{f_tm}C | R:{r_tm}C"))
            db.session.commit()
    return render_template('primer.html', primer_results=primer_results, input_seq=input_seq)

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
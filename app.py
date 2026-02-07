import os
from flask import Flask, render_template, request
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def index():
    results = None
    if request.method == "POST":
        # Get sequence and clean it (remove newlines, spaces, and numbers)
        raw_seq = request.form.get("sequence", "").upper()
        clean_seq = "".join(raw_seq.split())
        clean_seq = "".join([char for char in clean_seq if char in "ATCGU"])
        
        # Replace U with T for DNA processing if user pastes RNA
        clean_seq = clean_seq.replace("U", "T")

        if clean_seq:
            my_seq = Seq(clean_seq)
            
            # Calculations
            results = {
                "length": len(my_seq),
                "gc_content": round(gc_fraction(my_seq) * 100, 2),
                "mol_weight": round(my_seq.complement().transcribe().translate().molecular_weight(), 2) if len(my_seq) >= 3 else 0,
                "translation": str(my_seq.translate(to_stop=True)),
                "complement": str(my_seq.complement()),
                "orfs": [i for i in range(len(clean_seq)-2) if clean_seq[i:i+3] == "ATG"]
            }

    return render_template("index.html", results=results)

if __name__ == "__main__":
    # Render dynamic port assignment
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
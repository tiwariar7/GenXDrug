import streamlit as st
import requests
import time
import json
import os
from pathlib import Path
import py3Dmol
from stmol import showmol
import tempfile
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import QED
import safe as sf
import datamol as dm
from urllib3.util import Retry
from requests import Session
from requests.adapters import HTTPAdapter
import zipfile
import io
import glob
import plotly.graph_objects as go
import plotly.express as px
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import google.generativeai as genai
from google.generativeai.types import HarmCategory, HarmBlockThreshold
from config import API_KEYS

# GenMol Generator Class
class GenMol_Generator:
    __default_params__ = {
        "num_molecules": 10,
        "temperature": 1.0,
        "noise": 0.0,
        'step_size': 1,
        'unique': True,
        'scoring': 'QED'
    }

    def __init__(self, invoke_url = 'https://health.api.nvidia.com/v1/biology/nvidia/genmol/generate', auth = None, **kwargs):
        self.invoke_url = invoke_url
        self.auth = auth
        self.session = Session()
        self.num_generate = kwargs.get('num_generate', 1)
        self.verbose = False
        self.max_retries = kwargs.get('max_retries', 5)
        self.retries = Retry(
            total = self.max_retries,
            backoff_factor = 0.1,
            status_forcelist = [400],
            allowed_methods = {'POST'},
        )
        self.session.mount(self.invoke_url, HTTPAdapter(max_retries = self.retries))

    def produce(self, molecules, num_generate):
        generated = []

        for m in molecules:
            safe_segs = m.split('.')
            pos = np.random.randint(len(safe_segs))
            safe_segs[pos] = '[*{%d-%d}]' % (len(safe_segs[pos]), len(safe_segs[pos]) + 5)
            smiles = '.'.join(safe_segs)

            new_molecules = self.inference(
                smiles = smiles,
                num_molecules = max(10, num_generate),
                temperature = 1.5,
                noise = 2.0
            )

            new_molecules = [_['smiles'] for _ in new_molecules]

            if len(new_molecules) == 0:
                return []

            new_molecules = new_molecules[:(min(self.num_generate, len(new_molecules)))]
            generated.extend(new_molecules)

        self.molecules = list(set(generated))
        return self.molecules

    def inference(self, **params):
        headers = {
            "Authorization": "" if self.auth is None else "Bearer " + self.auth,
            "Content-Type": "application/json"
        }

        task = GenMol_Generator.__default_params__.copy()
        task.update(params)

        if self.verbose:
            print("TASK:", str(task))

        json_data = {k : str(v) for k, v in task.items()}

        response = self.session.post(self.invoke_url, headers=headers, json=json_data)
        response.raise_for_status()

        output = response.json()
        assert output['status'] == 'success'
        return output['molecules']

# Utils class for molecule manipulation
class Utils:
    @staticmethod
    def smiles2safe(smiles):
        converter = sf.SAFEConverter()
        return converter.encoder(smiles, allow_empty=True)

    @staticmethod
    def safe2smiles(safe):
        converter = sf.SAFEConverter()
        return converter.decoder(safe)

# Set page config
st.set_page_config(
    page_title="NVIDIA Bio AI Tools",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Configure Gemini
GEMINI_API_URL = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key={API_KEYS['gemini']}"

# Function to get AI response
def get_ai_response(user_query, context=""):
    try:
        prompt = f"""
        You are a helpful AI assistant for the NVIDIA Bio AI Tools application. 
        The application has the following tools:
        1. Protein Structure Prediction (AlphaFold2)
        2. Molecular Generation (GenMol)
        3. BioNeMo (Biological Sequence Analysis)
        4. Molecular Docking (DiffDock)
        
        Context: {context}
        User Query: {user_query}
        
        Please provide a helpful and accurate response. If the query is about a specific tool, 
        provide detailed information about that tool's functionality and usage.
        """
        
        # Prepare the request data
        data = {
            "contents": [{
                "parts": [{"text": prompt}]
            }]
        }
        
        # Make the API request
        headers = {
            'Content-Type': 'application/json'
        }
        
        response = requests.post(
            GEMINI_API_URL,
            headers=headers,
            json=data
        )
        
        # Check if the request was successful
        if response.status_code == 200:
            result = response.json()
            if 'candidates' in result and len(result['candidates']) > 0:
                if 'content' in result['candidates'][0]:
                    parts = result['candidates'][0]['content'].get('parts', [])
                    if parts:
                        return parts[0].get('text', 'No response generated')
            return "No response generated"
        else:
            return f"Error: {response.status_code} - {response.text}"
            
    except Exception as e:
        return f"Error getting AI response: {str(e)}"

# Sidebar navigation
st.sidebar.title("GenXDrug : Advancing Research with Next-Gen Drug Solution")
page = st.sidebar.radio(
    "Select Tool",
    ["Protein Structure Prediction", "Molecular Generation", "BioNeMo", "Molecular Docking", "AI Assistant"]
)

# Common functions
def clean_api_key(key):
    if key is None:
        return ""
    return re.sub(r'\s+', '', key)

# Protein Structure Prediction Page
if page == "Protein Structure Prediction":
    st.title("Protein Structure Prediction")
    st.markdown("""
    This tool allows you to predict protein structures using either AlphaFold2 or OpenFold. You can:
    - Input a protein sequence
    - Choose between AlphaFold2 and OpenFold
    - Predict its 3D structure
    - Visualize the predicted structure
    - Download the PDB file
    """)

    # Model selection
    model = st.radio(
        "Select Prediction Model",
        ["AlphaFold2", "OpenFold"],
        horizontal=True
    )

    # Sequence input with validation
    def validate_sequence(seq):
        seq = re.sub(r'\s+', '', seq)
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        if not all(aa in valid_aa for aa in seq):
            return False, "Sequence contains invalid amino acids"
        if len(seq) < 10:
            return False, "Sequence is too short (minimum 10 amino acids)"
        if len(seq) > 2000:
            return False, "Sequence is too long (maximum 2000 amino acids)"
        return True, seq

    sequence = st.text_area(
        "Protein Sequence",
        value="GGSKENEISHHAKEIERLQKEIERHKQSIKKLKQSEQSNPPPNPEGTRQARRNRRRRWRERQRQKENEISHHAKEIERLQKEIERHKQSIKKLKQSEC",
        height=200
    )

    # Advanced parameters
    with st.expander("Advanced Parameters"):
        if model == "AlphaFold2":
            algorithm = st.selectbox("Algorithm", ["mmseqs2"], index=0)
            e_value = st.number_input("E-value", value=0.0001, format="%.4f")
            iterations = st.number_input("Iterations", value=1, min_value=1, max_value=5)
            relax_prediction = st.checkbox("Relax Prediction", value=False)
            skip_template_search = st.checkbox("Skip Template Search", value=True)
        else:  # OpenFold
            st.markdown("### OpenFold Parameters")
            selected_models = st.multiselect(
                "Select Models",
                [1, 2, 3, 4, 5],
                default=[1, 2]
            )
            relax_prediction = st.checkbox("Relax Prediction", value=False)
            
            st.markdown("### MSA Alignment (Optional)")
            # Create MSA sequence that matches the input sequence length
            if sequence:
                seq_length = len(re.sub(r'\s+', '', sequence))
                msa_sequence = f">BQXYMDHSRWGGVPIWVK\n{sequence}\n>A0A076V4A1_9HIV1\n{'-' * seq_length}\n>BQXYMDHSRWGGVPIWVK\n{sequence}\n>A0A076V4A1_9HIV1\n{'-' * seq_length}\n>BQXYMDHSRWGGVPIWVK\n{sequence}\n>A0A076V4A1_9HIV1\n{'-' * seq_length}"
            else:
                msa_sequence = ">BQXYMDHSRWGGVPIWVK\nGGSKENEISHHAKEIERLQKEIERHKQSIKKLKQSEQSNPPPNPEGTRQARRNRRRRWRERQRQKENEISHHAKEIERLQKEIERHKQSIKKLKQSEC\n>A0A076V4A1_9HIV1\n----------------------------------------------------------------------\n>BQXYMDHSRWGGVPIWVK\nGGSKENEISHHAKEIERLQKEIERHKQSIKKLKQSEQSNPPPNPEGTRQARRNRRRRWRERQRQKENEISHHAKEIERLQKEIERHKQSIKKLKQSEC\n>A0A076V4A1_9HIV1\n----------------------------------------------------------------------\n>BQXYMDHSRWGGVPIWVK\nGGSKENEISHHAKEIERLQKEIERHKQSIKKLKQSEQSNPPPNPEGTRQARRNRRRRWRERQRQKENEISHHAKEIERLQKEIERHKQSIKKLKQSEC\n>A0A076V4A1_9HIV1\n----------------------------------------------------------------------"

            msa_sequence = st.text_area(
                "MSA Sequence",
                value=msa_sequence,
                height=150,
                help="Enter MSA alignment in A3M format. The sequence length must match the input protein sequence."
            )

            # Add validation for MSA sequence
            if model == "OpenFold" and msa_sequence:
                msa_lines = msa_sequence.split('\n')
                if len(msa_lines) < 2:
                    st.error("MSA sequence must contain at least one sequence")
                    st.stop()
                for i in range(0, len(msa_lines), 2):
                    if i + 1 >= len(msa_lines):
                        st.error("Invalid MSA format: missing sequence after header")
                        st.stop()
                    if not msa_lines[i].startswith('>'):
                        st.error("Invalid MSA format: sequence header must start with '>'")
                        st.stop()
                    if len(msa_lines[i + 1]) != seq_length:
                        st.error(f"MSA sequence length ({len(msa_lines[i + 1])}) does not match input sequence length ({seq_length})")
                        st.stop()

    if st.button("Predict Structure"):
        if not sequence:
            st.error("Please enter a protein sequence")
        else:
            is_valid, validation_result = validate_sequence(sequence)
            if not is_valid:
                st.error(validation_result)
                st.stop()
            
            with st.spinner("Predicting structure..."):
                try:
                    if model == "AlphaFold2":
                        url = "https://health.api.nvidia.com/v1/biology/deepmind/alphafold2"
                        data = {
                            "sequence": validation_result,
                            "algorithm": algorithm,
                            "e_value": e_value,
                            "iterations": iterations,
                            "databases": ["small_bfd"],
                            "relax_prediction": relax_prediction,
                            "skip_template_search": skip_template_search
                        }
                        headers = {
                            "content-type": "application/json",
                            "Authorization": "Bearer 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",
                            "NVCF-POLL-SECONDS": "300",
                        }
                    else:  # OpenFold
                        url = "https://health.api.nvidia.com/v1/biology/openfold/openfold2/predict-structure-from-msa-and-template"
                        data = {
                            "sequence": validation_result,
                            "alignments": {
                                "small_bfd": {
                                    "a3m": {
                                        "alignment": msa_sequence,
                                        "format": "a3m",
                                    }
                                },
                            },
                            "selected_models": selected_models,
                            "relax_prediction": relax_prediction,
                        }
                        headers = {
                            "content-type": "application/json",
                            "Authorization": "Bearer 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",
                            "NVCF-POLL-SECONDS": "300",
                        }

                    # Make the API request
                    response = requests.post(url, headers=headers, json=data, timeout=30)
                    
                    # Save the response to output1.json (but don't show it)
                    output_file = Path("output1.json")
                    output_file.write_text(response.text)
                    
                    if response.status_code == 200:
                        try:
                            response_data = response.json()
                            if model == "OpenFold":
                                if "structures_in_ranked_order" in response_data:
                                    pdb_string = response_data["structures_in_ranked_order"][0]["structure"]
                                else:
                                    st.error("Unexpected OpenFold response format")
                                    st.stop()
                            else:  # AlphaFold2
                                if isinstance(response_data, list) and len(response_data) > 0:
                                    pdb_string = response_data[0]
                                else:
                                    st.error("Unexpected AlphaFold2 response format")
                                    st.stop()
                        except json.JSONDecodeError as e:
                            st.error(f"Error parsing response: {str(e)}")
                            st.stop()
                    elif response.status_code == 202:
                        req_id = response.headers.get("nvcf-reqid")
                        if not req_id:
                            st.error("No request ID in response")
                            st.stop()
                            
                        status_url = "https://health.api.nvidia.com/v1/status"
                        max_attempts = 60
                        attempt = 0
                        while attempt < max_attempts:
                            status_response = requests.get(f"{status_url}/{req_id}", headers=headers, timeout=30)
                            if status_response.status_code != 202:
                                try:
                                    status_data = status_response.json()
                                    if isinstance(status_data, list) and len(status_data) > 0:
                                        pdb_string = status_data[0]
                                        break
                                    else:
                                        st.error("Unexpected status response format")
                                        st.stop()
                                except json.JSONDecodeError as e:
                                    st.error(f"Error parsing status response: {str(e)}")
                                    st.stop()
                            attempt += 1
                            time.sleep(5)
                        if attempt >= max_attempts:
                            st.error("Prediction timed out")
                            st.stop()
                    else:
                        st.error(f"Error: {response.status_code} - {response.text}")
                        st.stop()

                    # Create 3D visualization
                    st.subheader("Predicted Structure")
                    view = py3Dmol.view(width=1200, height=900)
                    view.addModel(pdb_string, "pdb")
                    view.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': 40, 'max': 100}}})
                    view.zoomTo()
                    showmol(view, height=900, width=1200)

                    # Download button
                    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
                        tmp.write(pdb_string.encode())
                        tmp_path = tmp.name

                    with open(tmp_path, "rb") as file:
                        st.download_button(
                            label="Download PDB File",
                            data=file,
                            file_name=f"predicted_structure_{model.lower()}.pdb",
                            mime="chemical/x-pdb"
                        )

                except Exception as e:
                    st.error(f"An error occurred: {str(e)}")
                    st.write("Full error details:", str(e.__class__.__name__))
                    import traceback
                    st.write("Traceback:", traceback.format_exc())

# Molecular Generation Page
elif page == "Molecular Generation":
    st.title("Molecular Generation")
    st.markdown("""
    This tool allows you to generate new molecules using NVIDIA's GenMol API. You can:
    - Generate molecules from scratch
    - Generate molecules based on a starting SMILES
    - Visualize the generated molecules
    - Analyze generation statistics
    """)

    # Use API key from config
    api_key = API_KEYS["genmol"]

    # Input options
    generation_type = st.radio(
        "Generation Type",
        ["From Scratch", "From SMILES"]
    )

    if generation_type == "From SMILES":
        input_smiles = st.text_input(
            "Input SMILES",
            value="CC(C=C(C(F)(F)F)C=C1OCC2=CC=CC=C2)=C1C3=CC=C(N(N=O)[C@@H]4CCCN(C)C4)N=N3",
            help="Enter a valid SMILES string"
        )
        
        # Validate SMILES
        try:
            testmol = dm.to_mol(input_smiles)
            if testmol is None:
                st.error("Invalid SMILES string. Please enter a valid SMILES.")
                st.stop()
        except Exception as e:
            st.error(f"Error validating SMILES: {str(e)}")
            st.stop()
    else:
        input_smiles = "[*{15-25}]"

    num_molecules = st.slider("Number of Molecules", 1, 50, 10)
    temperature = st.slider("Temperature", 0.1, 2.0, 1.0)
    noise = st.slider("Noise", 0.0, 2.0, 0.0)
    unique = st.checkbox("Generate Unique Molecules", True)

    if st.button("Generate Molecules"):
        if not api_key:
            st.error("Please enter your API key")
        elif not api_key.startswith("nvapi-"):
            st.error("API key must start with 'nvapi-'")
        else:
            with st.spinner("Generating molecules..."):
                try:
                    # Initialize GenMol generator
                    generator = GenMol_Generator(
                        auth=clean_api_key(api_key),
                        num_generate=num_molecules,
                        max_retries=5
                    )

                    # Generate molecules
                    if generation_type == "From SMILES":
                        # Convert SMILES to SAFE format
                        safe_string = Utils.smiles2safe(input_smiles)
                        # Generate new molecules
                        molecules = generator.produce([safe_string], num_molecules)
                    else:
                        # Generate from scratch
                        molecules = generator.inference(
                            smiles=input_smiles,
                            num_molecules=num_molecules,
                            temperature=temperature,
                            noise=noise
                        )

                    if not molecules:
                        st.error("No molecules were generated. Please try again with different parameters.")
                        st.stop()

                    # Extract SMILES strings from the response
                    if isinstance(molecules[0], dict):
                        smiles_list = [m['smiles'] for m in molecules]
                    else:
                        smiles_list = molecules

                    # Convert SAFE to SMILES if needed
                    if isinstance(smiles_list[0], str) and '[*{' in smiles_list[0]:
                        smiles_list = [sf.SAFEConverter(ignore_stereo=True).decoder(m) for m in smiles_list]

                    # Visualize molecules
                    ms = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
                    img = Draw.MolsToGridImage(
                        ms,
                        molsPerRow=4,
                        subImgSize=(240, 150),
                        legends=[f"Molecule {i+1}" for i in range(len(smiles_list))]
                    )
                    st.image(img)

                    # Calculate QED scores
                    qed_scores = [QED.qed(m) for m in ms]
                    
                    st.subheader("Generation Statistics")
                    df = pd.DataFrame({
                        'Molecule': [f"Molecule {i+1}" for i in range(len(smiles_list))],
                        'QED Score': qed_scores
                    })
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write("Score Statistics")
                        st.write(df.describe())
                    
                    with col2:
                        st.write("Score Distribution")
                        fig, ax = plt.subplots()
                        df['QED Score'].hist(ax=ax, bins=10)
                        st.pyplot(fig)

                    st.subheader("Generated SMILES")
                    for i, (smiles, score) in enumerate(zip(smiles_list, qed_scores)):
                        st.write(f"{i+1}. {smiles} (QED: {score:.3f})")

                except Exception as e:
                    st.error(f"An error occurred: {str(e)}")
                    st.write("Full error details:", str(e.__class__.__name__))
                    import traceback
                    st.write("Traceback:", traceback.format_exc())

# BioNeMo Page
elif page == "BioNeMo":
    st.title("BioNeMo: Biological Sequence Analysis")
    st.markdown("""
    This tool allows you to analyze biological sequences using NVIDIA's BioNeMo. You can:
    - Input DNA, RNA, or protein sequences
    - Perform sequence analysis
    - Get detailed sequence information
    - Analyze sequence motifs
    - Visualize structures and alignments
    """)

    # Use API key from config
    api_key = API_KEYS["bionemo"]

    # Sequence type selection
    sequence_type = st.selectbox(
        "Sequence Type",
        ["DNA", "RNA", "Protein"]
    )

    # Sequence input
    sequence = st.text_area(
        "Enter Sequence",
        height=200,
        help="Enter your biological sequence (DNA, RNA, or protein)"
    )

    # Enhanced Analysis Options
    with st.expander("Analysis Options"):
        analyze_gc_content = st.checkbox("GC Content Analysis", value=True)
        analyze_amino_acid = st.checkbox("Amino Acid Composition", value=True)
        analyze_secondary_structure = st.checkbox("Secondary Structure Prediction", value=True)
        analyze_motifs = st.checkbox("Motif Analysis", value=True)
        visualize_3d = st.checkbox("3D Structure Visualization", value=True)
        perform_alignment = st.checkbox("Sequence Alignment", value=True)
        analyze_phylogeny = st.checkbox("Phylogenetic Analysis", value=True)
        
        if analyze_motifs:
            motif_options = st.multiselect(
                "Select Motif Types to Analyze",
                ["DNA Binding Sites", "Protein Domains", "Regulatory Elements", "Conserved Regions"],
                default=["DNA Binding Sites", "Protein Domains"]
            )
            min_motif_length = st.slider("Minimum Motif Length", 3, 20, 5)
            max_motif_length = st.slider("Maximum Motif Length", 5, 50, 15)
            min_occurrences = st.slider("Minimum Occurrences", 1, 10, 2)

        if perform_alignment:
            st.write("### Multiple Sequence Input")
            num_sequences = st.number_input("Number of sequences to compare", min_value=2, max_value=10, value=2)
            sequences = []
            for i in range(num_sequences):
                seq = st.text_area(f"Sequence {i+1}", height=100)
                if seq:
                    sequences.append(seq)

    if st.button("Analyze Sequence"):
        if not api_key:
            st.error("Please enter your API key")
        elif not sequence:
            st.error("Please enter a sequence")
        else:
            with st.spinner("Analyzing sequence..."):
                try:
                    # Clean the sequence
                    sequence = re.sub(r'\s+', '', sequence).upper()
                    
                    # Basic sequence validation
                    if sequence_type == "DNA":
                        if not all(base in 'ATCG' for base in sequence):
                            st.error("Invalid DNA sequence. Only A, T, C, G are allowed.")
                            st.stop()
                    elif sequence_type == "RNA":
                        if not all(base in 'AUCG' for base in sequence):
                            st.error("Invalid RNA sequence. Only A, U, C, G are allowed.")
                            st.stop()
                    elif sequence_type == "Protein":
                        if not all(aa in 'ACDEFGHIKLMNPQRSTVWY' for aa in sequence):
                            st.error("Invalid protein sequence. Only standard amino acids are allowed.")
                            st.stop()

                    # Display basic sequence information
                    st.subheader("Sequence Information")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"Sequence Type: {sequence_type}")
                        st.write(f"Sequence Length: {len(sequence)}")
                        if sequence_type in ["DNA", "RNA"]:
                            st.write(f"Molecular Weight: {len(sequence) * 330:.2f} g/mol")
                        else:
                            st.write(f"Molecular Weight: {len(sequence) * 110:.2f} g/mol")

                    # GC Content Analysis
                    if analyze_gc_content and sequence_type in ["DNA", "RNA"]:
                        gc_count = sequence.count('G') + sequence.count('C')
                        gc_content = (gc_count / len(sequence)) * 100
                        with col2:
                            st.write(f"GC Content: {gc_content:.2f}%")
                            st.write(f"AT Content: {100 - gc_content:.2f}%")

                    # Amino Acid Composition
                    if analyze_amino_acid and sequence_type == "Protein":
                        aa_counts = {}
                        for aa in 'ACDEFGHIKLMNPQRSTVWY':
                            aa_counts[aa] = sequence.count(aa)
                        
                        st.subheader("Amino Acid Composition")
                        fig, ax = plt.subplots(figsize=(10, 4))
                        ax.bar(aa_counts.keys(), aa_counts.values())
                        ax.set_xlabel('Amino Acid')
                        ax.set_ylabel('Count')
                        ax.set_title('Amino Acid Distribution')
                        st.pyplot(fig)

                    # Secondary Structure Prediction
                    if analyze_secondary_structure and sequence_type == "Protein":
                        st.subheader("Secondary Structure Prediction")
                        # This is a placeholder for actual secondary structure prediction
                        # In a real implementation, you would call the BioNeMo API here
                        st.info("Secondary structure prediction would be implemented here using BioNeMo API")

                    # Motif Analysis
                    if analyze_motifs:
                        st.subheader("Motif Analysis")
                        
                        # Prepare data for motif analysis
                        headers = {
                            "Authorization": f"Bearer {api_key}",
                            "Content-Type": "application/json"
                        }
                        
                        data = {
                            "sequence": sequence,
                            "sequence_type": sequence_type,
                            "motif_types": motif_options,
                            "min_length": min_motif_length,
                            "max_length": max_motif_length,
                            "min_occurrences": min_occurrences
                        }
                        
                        try:
                            # Call BioNeMo API for motif analysis
                            response = requests.post(
                                "https://health.api.nvidia.com/v1/biology/bionemo/motif-analysis",
                                headers=headers,
                                json=data
                            )
                            response.raise_for_status()
                            motif_results = response.json()
                            
                            # Display motif analysis results
                            for motif_type in motif_options:
                                st.write(f"### {motif_type}")
                                
                                if motif_type in motif_results:
                                    motifs = motif_results[motif_type]
                                    if motifs:
                                        # Create a DataFrame for better visualization
                                        motif_data = []
                                        for motif in motifs:
                                            motif_data.append({
                                                "Motif": motif["sequence"],
                                                "Position": f"{motif['start']}-{motif['end']}",
                                                "Score": f"{motif['score']:.2f}",
                                                "Significance": motif["significance"]
                                            })
                                        
                                        df = pd.DataFrame(motif_data)
                                        st.dataframe(df)
                                        
                                        # Visualize motif positions
                                        fig, ax = plt.subplots(figsize=(10, 2))
                                        for motif in motifs:
                                            ax.plot([motif["start"], motif["end"]], [0, 0], 
                                                   linewidth=10, solid_capstyle='butt')
                                        ax.set_xlim(0, len(sequence))
                                        ax.set_ylim(-0.5, 0.5)
                                        ax.set_yticks([])
                                        ax.set_xlabel("Sequence Position")
                                        ax.set_title(f"{motif_type} Positions")
                                        st.pyplot(fig)
                                    else:
                                        st.info(f"No significant {motif_type} found in the sequence.")
                                else:
                                    st.warning(f"Analysis for {motif_type} is not available.")
                            
                        except requests.exceptions.RequestException as e:
                            st.error(f"Error in motif analysis: {str(e)}")
                            st.info("Using placeholder data for demonstration")
                            
                            # Placeholder data for demonstration
                            placeholder_motifs = {
                                "DNA Binding Sites": [
                                    {"sequence": "TATAAA", "start": 50, "end": 55, "score": 0.85, "significance": "High"},
                                    {"sequence": "GCGCGC", "start": 120, "end": 125, "score": 0.75, "significance": "Medium"}
                                ],
                                "Protein Domains": [
                                    {"sequence": "LHLH", "start": 30, "end": 33, "score": 0.90, "significance": "High"},
                                    {"sequence": "EFH", "start": 80, "end": 82, "score": 0.80, "significance": "Medium"}
                                ]
                            }
                            
                            for motif_type in motif_options:
                                if motif_type in placeholder_motifs:
                                    st.write(f"### {motif_type} (Placeholder Data)")
                                    motifs = placeholder_motifs[motif_type]
                                    df = pd.DataFrame(motifs)
                                    st.dataframe(df)
                                    
                                    # Visualize motif positions
                                    fig, ax = plt.subplots(figsize=(10, 2))
                                    for motif in motifs:
                                        ax.plot([motif["start"], motif["end"]], [0, 0], 
                                               linewidth=10, solid_capstyle='butt')
                                    ax.set_xlim(0, len(sequence))
                                    ax.set_ylim(-0.5, 0.5)
                                    ax.set_yticks([])
                                    ax.set_xlabel("Sequence Position")
                                    ax.set_title(f"{motif_type} Positions (Placeholder)")
                                    st.pyplot(fig)

                    # Display sequence with formatting
                    st.subheader("Formatted Sequence")
                    if sequence_type == "DNA":
                        formatted_seq = ' '.join([sequence[i:i+10] for i in range(0, len(sequence), 10)])
                        st.code(formatted_seq)
                    elif sequence_type == "RNA":
                        formatted_seq = ' '.join([sequence[i:i+10] for i in range(0, len(sequence), 10)])
                        st.code(formatted_seq)
                    else:  # Protein
                        formatted_seq = ' '.join([sequence[i:i+10] for i in range(0, len(sequence), 10)])
                        st.code(formatted_seq)

                    # Enhanced Visualization Section
                    if visualize_3d and sequence_type == "Protein":
                        st.subheader("3D Structure Visualization")
                        try:
                            # Call BioNeMo API for structure prediction
                            headers = {
                                "Authorization": f"Bearer {api_key}",
                                "Content-Type": "application/json"
                            }
                            
                            data = {
                                "sequence": sequence,
                                "model": "alphafold2"  # or other structure prediction model
                            }
                            
                            response = requests.post(
                                "https://health.api.nvidia.com/v1/biology/bionemo/structure-prediction",
                                headers=headers,
                                json=data
                            )
                            response.raise_for_status()
                            structure_data = response.json()
                            
                            # Create 3D visualization using py3Dmol
                            view = py3Dmol.view(width=800, height=600)
                            view.addModel(structure_data["pdb"], "pdb")
                            view.setStyle({'model': 0}, {'cartoon': {'color': 'spectrum'}})
                            view.addSurface(py3Dmol.VDW, {'opacity': 0.4, 'color': 'white'})
                            view.zoomTo()
                            showmol(view, height=600, width=800)
                            
                            # Add structure analysis
                            st.write("### Structure Analysis")
                            st.write(f"Confidence Score: {structure_data.get('confidence_score', 'N/A')}")
                            st.write(f"Predicted Domains: {structure_data.get('domains', 'N/A')}")
                            
                        except Exception as e:
                            st.error(f"Error in structure prediction: {str(e)}")
                            st.info("Using placeholder structure for demonstration")
                            # Placeholder visualization code here

                    # Sequence Alignment Visualization
                    if perform_alignment and len(sequences) >= 2:
                        st.subheader("Sequence Alignment")
                        try:
                            # Create multiple sequence alignment
                            alignment = MultipleSeqAlignment([
                                SeqRecord(Seq(seq), id=f"Seq{i+1}") 
                                for i, seq in enumerate(sequences)
                            ])
                            
                            # Create alignment visualization
                            fig = go.Figure()
                            
                            # Add sequence letters as heatmap
                            alignment_array = np.array([list(str(rec.seq)) for rec in alignment])
                            fig.add_trace(go.Heatmap(
                                z=alignment_array,
                                colorscale='Viridis',
                                showscale=False
                            ))
                            
                            # Customize layout
                            fig.update_layout(
                                title="Multiple Sequence Alignment",
                                xaxis_title="Position",
                                yaxis_title="Sequence",
                                height=400
                            )
                            
                            st.plotly_chart(fig, use_container_width=True)
                            
                            # Add conservation analysis
                            conservation = []
                            for i in range(len(alignment[0])):
                                column = alignment_array[:, i]
                                unique = np.unique(column)
                                conservation.append(len(unique) == 1)
                            
                            # Plot conservation
                            fig_cons = go.Figure()
                            fig_cons.add_trace(go.Scatter(
                                y=conservation,
                                mode='lines',
                                name='Conservation'
                            ))
                            fig_cons.update_layout(
                                title="Sequence Conservation",
                                xaxis_title="Position",
                                yaxis_title="Conserved",
                                height=200
                            )
                            st.plotly_chart(fig_cons, use_container_width=True)
                            
                        except Exception as e:
                            st.error(f"Error in sequence alignment: {str(e)}")

                    # Phylogenetic Analysis
                    if analyze_phylogeny and len(sequences) >= 2:
                        st.subheader("Phylogenetic Analysis")
                        try:
                            # Create distance matrix
                            calculator = DistanceCalculator('identity')
                            constructor = DistanceTreeConstructor(calculator)
                            
                            # Create alignment
                            alignment = MultipleSeqAlignment([
                                SeqRecord(Seq(seq), id=f"Seq{i+1}") 
                                for i, seq in enumerate(sequences)
                            ])
                            
                            # Build tree
                            tree = constructor.build_tree(alignment)
                            
                            # Create phylogenetic tree visualization
                            fig = go.Figure()
                            
                            # Add tree branches
                            for clade in tree.get_nonterminals():
                                for child in clade:
                                    fig.add_trace(go.Scatter(
                                        x=[clade.branch_length, child.branch_length],
                                        y=[clade.name, child.name],
                                        mode='lines',
                                        line=dict(width=2)
                                    ))
                            
                            # Customize layout
                            fig.update_layout(
                                title="Phylogenetic Tree",
                                showlegend=False,
                                height=600
                            )
                            
                            st.plotly_chart(fig, use_container_width=True)
                            
                            # Add distance matrix heatmap
                            distances = calculator.get_distance(alignment)
                            fig_dist = go.Figure(data=go.Heatmap(
                                z=distances,
                                x=[f"Seq{i+1}" for i in range(len(sequences))],
                                y=[f"Seq{i+1}" for i in range(len(sequences))],
                                colorscale='Viridis'
                            ))
                            
                            fig_dist.update_layout(
                                title="Distance Matrix",
                                height=400
                            )
                            
                            st.plotly_chart(fig_dist, use_container_width=True)
                            
                        except Exception as e:
                            st.error(f"Error in phylogenetic analysis: {str(e)}")

                except Exception as e:
                    st.error(f"An error occurred: {str(e)}")

# Add Molecular Docking Page
elif page == "Molecular Docking":
    st.title("Molecular Docking with DiffDock")
    st.markdown("""
    This tool allows you to perform molecular docking using NVIDIA's DiffDock. You can:
    - Upload protein (PDB) and ligand (SDF) files
    - Run docking simulations
    - Visualize docking results
    - Download docking poses
    """)

    # Use API key from config
    api_key = API_KEYS["diffdock"]

    # Helper functions for docking
    def file_to_json_compatible_string(file_content):
        return file_content.decode("utf-8")

    def run_diffdock(protein_content, ligand_content):
        data = {
            "ligand": file_to_json_compatible_string(ligand_content),
            "ligand_file_type": "sdf",
            "protein": file_to_json_compatible_string(protein_content),
            "num_poses": 10,
            "time_divisions": 20,
            "steps": 18,
            "save_trajectory": False,
            "is_staged": False
        }

        headers = {
            "Authorization": f'Bearer {api_key}',
            "Content-Type": "application/json",
            "accept": "application/json",
        }

        response = requests.post("https://health.api.nvidia.com/v1/biology/mit/diffdock", 
                               headers=headers, json=data)
        response.raise_for_status()
        return response.json()

    def create_3d_view(protein_content, pose_mol, pose_index, confidence_scores):
        view = py3Dmol.view(width=1200, height=900)
        
        # Add protein with enhanced styling
        view.addModel(protein_content.decode(), 'pdb')
        view.setStyle({'model': 0}, {'cartoon': {'color': 'white', 'opacity': 0.7}})
        view.setViewStyle({'style': 'outline', 'color': 'black', 'width': 0.03})
        view.addSurface(py3Dmol.VDW, {'opacity': 0.4, 'color': 'white'})

        # Add ligand with surface
        pose_block = Chem.MolToMolBlock(pose_mol)
        view.addModel(pose_block, 'mol')
        view.setStyle({'model': 1}, {'stick': {'radius': 0.3, 'colorscheme': 'magentaCarbon'}})
        view.addSurface(py3Dmol.VDW, {'opacity': 0.7, 'colorscheme': 'magentaCarbon'}, {'model': 1})

        # Zoom to complex
        view.zoomTo()
        
        # Get confidence score
        score = round(confidence_scores[pose_index], 3)
        score_color = "green" if score > -0.5 else "blue" if score >= -1.5 else "red"
        
        return view, score, score_color

    # File upload
    protein_file = st.file_uploader("Upload Protein PDB file", type=["pdb"])
    ligand_files = st.file_uploader("Upload Ligand SDF or ZIP of SDFs", type=["sdf", "zip"])

    if protein_file and ligand_files:
        with st.spinner("Processing files..."):
            # Handle ZIP upload
            if ligand_files.name.endswith(".zip"):
                with zipfile.ZipFile(ligand_files) as z:
                    sdf_files = [f for f in z.namelist() if f.endswith(".sdf")]
                    sdf_files.sort(key=lambda x: int(x.split("_")[1].split(".")[0]))
                    ligand_contents = {name: z.read(name) for name in sdf_files}
            else:
                ligand_contents = {ligand_files.name: ligand_files.read()}

            # Store in session state
            st.session_state.protein_content = protein_file.read()
            st.session_state.ligand_contents = ligand_contents

    if "ligand_contents" in st.session_state:
        ligand_selector = st.selectbox("Select Ligand", list(st.session_state.ligand_contents.keys()))
        
        if st.button("Run Docking"):
            with st.spinner("Running DiffDock..."):
                try:
                    result = run_diffdock(
                        st.session_state.protein_content,
                        st.session_state.ligand_contents[ligand_selector]
                    )
                    
                    # Store results in session state
                    st.session_state.results = {
                        "poses": [Chem.MolFromMolBlock(geom) for geom in result["ligand_positions"]],
                        "scores": result["position_confidence"]
                    }
                    
                    st.success("Docking completed!")
                except Exception as e:
                    st.error(f"Error running docking: {str(e)}")

    # Visualization and results
    if "results" in st.session_state:
        st.header("Docking Results")
        
        # Pose selection
        pose_idx = st.slider("Select Pose", 0, len(st.session_state.results["poses"])-1, 0)
        
        # Get current pose and score
        current_pose = st.session_state.results["poses"][pose_idx]
        current_score = st.session_state.results["scores"][pose_idx]
        
        # Create columns for layout
        col1, col2 = st.columns([3, 1])
        
        with col1:
            st.subheader("3D Visualization")
            
            # Create py3Dmol viewer
            view = py3Dmol.view(width=800, height=600)
            
            # Add protein with enhanced styling
            view.addModel(st.session_state.protein_content.decode(), "pdb")
            view.setStyle({'model': 0}, {'cartoon': {'color': 'white', 'opacity': 0.7}})
            view.setViewStyle({'style':'outline','color':'black','width':0.03})
            
            # Add protein surface
            view.addSurface(py3Dmol.VDW, {'opacity': 0.4, 'color': 'white'})
            
            # Add ligand with styling
            pose_block = Chem.MolToMolBlock(current_pose)
            view.addModel(pose_block, "mol")
            view.setStyle({'model': 1}, {'stick': {
                'radius': 0.3,
                'colorscheme': 'magentaCarbon'
            }})
            
            # Add ligand surface
            view.addSurface(py3Dmol.VDW, {
                'opacity': 0.7,
                'colorscheme': 'magentaCarbon'
            }, {'model': 1})
            
            view.zoomTo()
            showmol(view, height=600, width=800)

        with col2:
            st.subheader("Analysis")
            
            # Score display with color coding
            score_color = "green" if current_score > -0.5 else "orange" if current_score > -1.5 else "red"
            st.markdown(f"""
            **Ligand Name**: `{ligand_selector}`  
            **Confidence Score**: <span style='color:{score_color}; font-size: 24px'>{current_score:.2f}</span>
            """, unsafe_allow_html=True)
            
            # Pose statistics
            st.write("### Pose Statistics")
            st.write(f"Total Poses: {len(st.session_state.results['poses'])}")
            st.write(f"Current Pose: {pose_idx + 1}")
            
            # Download current pose
            pose_bytes = Chem.MolToMolBlock(current_pose).encode()
            st.download_button(
                label="Download Current Pose",
                data=pose_bytes,
                file_name=f"pose_{pose_idx+1}.sdf",
                mime="chemical/x-sdf"
            )
            
            # Download all poses
            with tempfile.TemporaryDirectory() as tmpdir:
                output_file = os.path.join(tmpdir, "results.zip")
                with zipfile.ZipFile(output_file, "w") as z:
                    for i, (pose, score) in enumerate(zip(st.session_state.results["poses"], st.session_state.results["scores"])):
                        pose_sdf = Chem.MolToMolBlock(pose)
                        pose_path = os.path.join(tmpdir, f"pose_{i}.sdf")
                        with open(pose_path, "w") as f:
                            f.write(pose_sdf)
                        z.write(pose_path, f"pose_{i}.sdf")
                
                with open(output_file, "rb") as f:
                    st.download_button(
                        label="Download All Poses",
                        data=f.read(),
                        file_name="docking_results.zip",
                        mime="application/zip"
                    )

# Add AI Assistant page
elif page == "AI Assistant":
    st.title("AI Assistant")
    st.markdown("""
    Welcome to the AI Assistant! I can help you with:
    - Understanding the tools and their functionalities
    - Providing guidance on using specific features
    - Answering questions about bioinformatics
    - Troubleshooting issues
    """)
    
    # Initialize chat history in session state
    if "messages" not in st.session_state:
        st.session_state.messages = []
    
    # Display chat messages
    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])
    
    # Chat input
    if prompt := st.chat_input("Ask me anything about the tools or bioinformatics"):
        # Add user message to chat history
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)
        
        # Get context based on current tool
        context = ""
        if "current_tool" in st.session_state:
            context = f"Current tool: {st.session_state.current_tool}"
        
        # Get AI response
        with st.chat_message("assistant"):
            with st.spinner("Thinking..."):
                response = get_ai_response(prompt, context)
                st.markdown(response)
                st.session_state.messages.append({"role": "assistant", "content": response})

# Update help information
st.sidebar.markdown("""
### Help
- **Protein Structure Prediction**: Predict 3D structures using AlphaFold2
- **Molecular Generation**: Generate new molecules using GenMol
- **BioNeMo**: Biological sequence analysis tools
- **Molecular Docking**: Perform protein-ligand docking using DiffDock
- **AI Assistant**: Get help and guidance from our AI assistant
""")

# Add current tool tracking
if page != "AI Assistant":
    st.session_state.current_tool = page 

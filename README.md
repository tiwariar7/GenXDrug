# GenXDrug - NVIDIA Bio AI Tools

![GenXDrug Logo](https://github.com/tiwariar7/GenXDrug/blob/main/Supporting%20Artifact/Architecture_Diagram-GenXDrug.svg)

## Overview

GenXDrug is a comprehensive web application that integrates multiple NVIDIA Bio AI tools for molecular and biological sequence analysis. The application provides a user-friendly interface for protein structure prediction, molecular generation, biological sequence analysis, and molecular docking.

## Live Demo

Check out the live demo of GenXDrug: [GenXDrug MVP](https://tiwariar7-genxdrug-app-afb80f.streamlit.app/)

## Features

### 1. Protein Structure Prediction (AlphaFold2)
- Input protein sequences for 3D structure prediction
- Visualize predicted protein structures
- Download PDB files
- Advanced parameter customization

### 2. Molecular Generation (GenMol)
- Generate novel molecules from scratch
- Generate molecules based on input SMILES
- Visualize generated molecules
- Analyze generation statistics
- Download generated structures

### 3. BioNeMo (Biological Sequence Analysis)
- Analyze DNA, RNA, and protein sequences
- Perform GC content analysis
- Analyze amino acid composition
- Predict secondary structures
- Identify sequence motifs
- Visualize sequence alignments
- Perform phylogenetic analysis

### 4. Molecular Docking (DiffDock)
- Upload protein and ligand files
- Run docking simulations
- Visualize docking results
- Analyze binding poses
- Download docking results

### 5. AI Assistant
- Get help and guidance
- Answer questions about tools
- Provide bioinformatics insights
- Troubleshoot issues

## Architecture

![Architecture Diagram](https://github.com/tiwariar7/GenXDrug/blob/main/Supporting%20Artifact/Architecture_Diagram-GenXDrug.svg)

## UI/UX Design

![Prototype Wireframing](https://github.com/tiwariar7/GenXDrug/blob/main/Supporting%20Artifact/Prototype_Wireframing-GenXDrug.svg)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/tiwariar7/genxdrug.git
cd genxdrug
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Create a `config.py` file with your API keys:
```python
API_KEYS = {
    "alphafold": "your_alphafold_api_key",
    "genmol": "your_genmol_api_key",
    "bionemo": "your_bionemo_api_key",
    "diffdock": "your_diffdock_api_key",
    "gemini": "your_gemini_api_key"
}
```

4. Run the application:
```bash
streamlit run app.py
```

## Troubleshooting

### Common Issues and Solutions

1. **API Key Errors**
   - **Issue**: "Invalid API key" or "API key not found" errors
   - **Solution**: 
     - Verify that your API keys are correctly formatted in `config.py`
     - Ensure there are no extra spaces or special characters
     - Check if the API keys are still valid and not expired
     - Make sure the `config.py` file is in the correct directory

2. **Dependency Installation Issues**
   - **Issue**: Error during `pip install -r requirements.txt`
   - **Solution**:
     - Try installing dependencies one by one
     - Use a virtual environment
     - For RDKit installation issues:
       ```bash
       conda install -c conda-forge rdkit
       ```
     - For Py3DMol issues:
       ```bash
       pip install py3dmol
       ```

3. **Streamlit Connection Issues**
   - **Issue**: "Connection refused" or "Cannot connect to Streamlit"
   - **Solution**:
     - Check if port 8501 is available
     - Try running with a different port:
       ```bash
       streamlit run app.py --server.port 8502
       ```
     - Ensure no firewall is blocking the connection

4. **Visualization Problems**
   - **Issue**: 3D molecule viewer not working
   - **Solution**:
     - Clear browser cache
     - Try a different browser
     - Ensure WebGL is enabled in your browser
     - Check if all required JavaScript libraries are loaded

5. **File Upload Issues**
   - **Issue**: Cannot upload PDB or SDF files
   - **Solution**:
     - Verify file format is correct
     - Check file size (should be under Streamlit's limit)
     - Ensure file is not corrupted
     - Try converting file to a different format if needed

6. **Memory Issues**
   - **Issue**: Application crashes or runs slowly
   - **Solution**:
     - Increase Streamlit's memory limit:
       ```bash
       streamlit run app.py --server.maxUploadSize=500
       ```
     - Close other memory-intensive applications
     - Use smaller input files when possible

7. **Browser Compatibility**
   - **Issue**: UI elements not displaying correctly
   - **Solution**:
     - Use the latest version of Chrome, Firefox, or Edge
     - Disable browser extensions that might interfere
     - Clear browser cache and cookies

8. **API Rate Limiting**
   - **Issue**: "Too many requests" error
   - **Solution**:
     - Wait before making new requests
     - Contact NVIDIA support for higher rate limits
     - Implement request queuing in your code

If you encounter any other issues, please:
1. Check the Streamlit logs for error messages
2. Verify your system meets the minimum requirements
3. Try restarting the application
4. Open an issue in the GitHub repository with:
   - Error message
   - Steps to reproduce
   - System information
   - Screenshots if applicable

## API Keys

The application requires the following API keys:
- NVIDIA AlphaFold2 API key
- NVIDIA GenMol API key
- NVIDIA BioNeMo API key
- NVIDIA DiffDock API key
- Google Gemini API key

## Dependencies

- streamlit
- requests
- numpy
- pandas
- matplotlib
- rdkit
- py3Dmol
- biopython
- plotly
- google-generativeai

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For any questions or suggestions, please open an issue in the GitHub repository.

## Acknowledgments

- NVIDIA for providing the Bio AI tools
- Streamlit for the web framework
- All contributors and users of the application 

import re
import requests
from bs4 import BeautifulSoup
from rdkit import Chem

# Fetch the web page
url = "https://arrows.emsl.pnnl.gov/api/queue_html/"  # Replace with the actual URL
response = requests.get(url)

# Print the response content for debugging
# print(response.content.decode('utf-8'))
# smiles_pattern = re.compile(r'([A-Za-z0-9@+\-\[\]\(\)=#$%.,]+)')

# Parse the web page content
soup = BeautifulSoup(response.content, "html.parser")

# Extract the preformatted text block
preformatted_block = soup.find("pre").text

# Split the block into lines
lines = preformatted_block.splitlines()

# Find the starting point of the table
start_index = None
for i, line in enumerate(lines):
    if "#chemdb_queue" in line:
        start_index = i
        break

if start_index is not None:
    # Skip the header lines and parse the data
    data = []
    for line in lines[
        start_index + 4 :
    ]:  # Skip the first four lines after the start point
        if line.strip():  # Skip empty lines
            # parts = line.split()
            parts = re.split(r"\s{2,}", line.strip(), maxsplit=3)
            if len(parts) == 4:  # Ensure there are enough parts
                queue_entry = parts[0]
                link = parts[1]
                machine_type = parts[2]
                esmiles = parts[3].strip()

                data.append(
                    {
                        "queue_entry": queue_entry,
                        "link": link,
                        "machine_type": machine_type,
                        # "fetched": fetched,
                        "esmiles": esmiles,
                        # "smiles": match.group(1) if match else None,
                    }
                )
            else:
                print(f"Skipping line due to unexpected format: {line}")


def extract_valid_smiles(input_string):
    # Define a regex pattern for candidate SMILES strings
    smiles_pattern = r"\b[a-zA-Z0-9@+\-\[\]\(\)=#$%]+\b"

    # Find all candidate substrings in the input string
    candidate_smiles = re.findall(smiles_pattern, input_string)

    # Validate each candidate using RDKit
    valid_smiles = [
        smiles for smiles in candidate_smiles if Chem.MolFromSmiles(smiles) is not None
    ]

    return valid_smiles


# Print the extracted data
for entry in data:
    if "censo" in entry["esmiles"].lower():
        print(entry["queue_entry"], extract_valid_smiles(entry["esmiles"]).pop())

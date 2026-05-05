import json 
from pathlib import Path 

default_decisions_path = "data/decisions.json"

def save_decisions(decisions, filepath = default_decisions_path): 
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)  
    with open(filepath, 'w') as f:   
        json.dump(decisions, f, indent = 2)



def load_decisions(filepath = default_decisions_path):
    filepath = Path(filepath)
    try:
        with open(filepath, 'r') as f:
            decisions = json.load(f)
            return decisions
    except FileNotFoundError:
        return []

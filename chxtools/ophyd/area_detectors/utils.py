import numpy as np

# Define a mapping between the name of a detector and a folder to which to write the data
# Needed for back-compatibility with old scripts using detector names different from the folder names
name_dir_mapping = {
        'xray_eye1': 'xray_eye-1',
        'xray_eye2': 'xray_eye-2',
        'xray_eye3': 'xray_eye-3',
        'xray_eye4': 'xray_eye-4',
        'OAV': 'oavcam-1',
        'BCam' : 'zoomcam-1',
        'eiger500k': 'eiger500k-1',
        'eiger500K_single': 'eiger500k-1',
        'eiger500k_single': 'eiger500k-1',
        'eiger1m': 'eiger1m-1',
        'eiger1m_single': 'eiger1m-1',
        'eiger4m': 'eiger4m-1',
        'eiger4m_single': 'eiger4m-1',
        'pilatus800': 'pilatus800k-1',
    }


def str_to_wav(string:str):
    return np.frombuffer(f"{string}\x00".encode('utf-8'),dtype="uint8")


def _asset_undefined_error():
    raise ValueError("File provider not specified, please initialize device with a path provider")
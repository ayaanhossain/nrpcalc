import atexit
import os
import shutil

current_uuid = None

def setup_proj_dir(uuid):
    global current_uuid
    current_uuid = uuid
    directory = './{}/'.format(current_uuid)
    if not os.path.exists(directory):
        os.makedirs(directory)

@atexit.register
def remove_proj_dir():
    global current_uuid
    directory = './{}/'.format(current_uuid)
    if not current_uuid is None and os.path.isdir(directory):
        shutil.rmtree(directory)
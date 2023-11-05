import sys

# import yaml

print(sys.executable)
print(sys.path)


import os

folder_path = r'D:\\mesh_data\\4\\outt'  # f g v

for filename in os.listdir(folder_path):
    # if filename.startswith('camel-') and filename.endswith('.obj'):
    if filename.startswith('out_'):
        new_name = filename.replace('out_', 'outt_')
        old_path = os.path.join(folder_path, filename)
        new_path = os.path.join(folder_path, new_name)
        os.rename(old_path, new_path)

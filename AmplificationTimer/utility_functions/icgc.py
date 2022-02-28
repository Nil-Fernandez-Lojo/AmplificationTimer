def icgc(samplename,config):
    file_white = config['path_folder_snv_white']/(samplename+config['suffix_snvs'])
    file_gray = config['path_folder_snv_gray']/(samplename+config['suffix_snvs'])
    return file_gray.exists() or file_white.exists()
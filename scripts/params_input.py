def parse_input(file_path):
    params = {}
    with open(file_path, 'r') as f:
        for line in f:
            # 去除前后空白字符
            line = line.strip()
            # 忽略注释行
            if not line.startswith('#'):
                # 分割键值对
                if '=' in line:
                    key, value = line.split('=', 1)
                    # 去除键和值的前后空白字符，并保存到字典中
                    params[key.strip()] = value.strip()
    return params

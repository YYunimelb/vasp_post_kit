import configparser

config = configparser.ConfigParser()

config.read("yy_input")
for i in (config.items("vasp")):
    print(i)
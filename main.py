import argparse
from scripts.params_input import parse_input

def main(config = {'GG': '11', 'TT': '11', 'HH': '11', 'YY': '22'}):


    print(config)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Your program description')

    parser.add_argument("--config_path", default="yy_input")
    envir_config = parse_input("config/configurations")



    args = parser.parse_args()
    config = parse_input(args.config_path)

    main(config,envir_config)

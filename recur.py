import sys
import os

recur_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'src'))
sys.path.insert(0, recur_dir)

if __name__ == "__main__":
    from recur.run.__main__ import main
    args = sys.argv[1:]
    main(args)

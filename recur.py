# import sys
# import os

# # Adjust the path to include the recur directory
# recur_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'src'))
# sys.path.insert(0, recur_dir)

# if __name__ == "__main__":
#     from recur.run.__main__ import main, setup_environment
#     if not os.getenv('ENV_SETUP_DONE'):
#         setup_environment()
#         os.environ['ENV_SETUP_DONE'] = '1'
#     args = sys.argv[1:]
#     main(args)


import sys
import os

# Add the src path to sys.path
base_path = os.path.abspath(os.path.dirname(__file__))
src_path = os.path.join(base_path, 'src')
sys.path.insert(0, src_path)

if __name__ == "__main__":
    from recur.run.__main__ import main, setup_environment
    if not os.getenv('ENV_SETUP_DONE'):
        setup_environment()
        os.environ['ENV_SETUP_DONE'] = '1'
    args = sys.argv[1:]
    main(args)

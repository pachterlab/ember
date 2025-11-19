from .cli import create_parser, run_command

def main():
    """Main entry point for the ember command-line tool."""
    parser = create_parser()
    args = parser.parse_args()
    run_command(parser, args)

if __name__ == "__main__":
    main()
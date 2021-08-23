from accuracode.accura.__init__ import __ASSAY__
from accuracode.tools.multi import Multi


class Multi_accura(Multi):
    pass


def main():
    multi = Multi_accura(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()

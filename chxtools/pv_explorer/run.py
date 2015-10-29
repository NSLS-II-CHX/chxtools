from __future__ import (unicode_literals, print_function, absolute_import, \
                        division)

from enaml.qt.qt_application import QtApplication
import enaml

from chxtools.pv_explorer.model import Model
with enaml.imports():
    from chxtools.pv_explorer.view import MainView

def main():
    app = QtApplication()
    model = Model()
    view = MainView(model=model)
    view.show()

    app.start()

if __name__ == "__main__":
    main()

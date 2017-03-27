from PyQt4 import QtCore, QtGui

class Widget(QtGui.QWidget):

    def __init__(self, *args, **kwargs):
        super(Widget, self).__init__(*args, **kwargs)
        self.resize(800,600)

        self.vlayout = QtGui.QVBoxLayout(self)
        self.table = QtGui.QTableView()
        self.vlayout.addWidget(self.table)

        self.hlayout = QtGui.QHBoxLayout()
        self.list1 = QtGui.QListView()
        self.list2 = QtGui.QListView()
        self.list3 = QtGui.QListView()
        self.hlayout.addWidget(self.list1)
        self.hlayout.addWidget(self.list2)
        self.hlayout.addWidget(self.list3)

        self.vlayout.addLayout(self.hlayout)

        self.model = QtGui.QStandardItemModel(10,10,self)
        self.table.setModel(self.model)

        self.list1.setModel(self.model)
        self.list1.setModelColumn(0)
        self.list2.setModel(self.model)
        self.list2.setModelColumn(1)
        self.list3.setModel(self.model)
        self.list3.setModelColumn(2)

        self.populateTable()

    def populateTable(self):
        for row in xrange(10):
            for col in xrange(10):
                item = QtGui.QStandardItem('%d-%d' % (row, col))
                self.model.setItem(row, col, item)


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication([])
    window = Widget()
    window.show()
    window.raise_()
    sys.exit(app.exec_())

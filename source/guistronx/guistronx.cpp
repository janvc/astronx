#include "guistronx_window.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    guistronx_window w;
    w.show();

    return a.exec();
}

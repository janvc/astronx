#ifndef GUISTRONX_WINDOW_H
#define GUISTRONX_WINDOW_H

#include <QMainWindow>

namespace Ui {
class guistronx_window;
}

class guistronx_window : public QMainWindow
{
    Q_OBJECT

public:
    explicit guistronx_window(QWidget *parent = 0);
    ~guistronx_window();

private:
    Ui::guistronx_window *ui;
};

#endif // GUISTRONX_WINDOW_H

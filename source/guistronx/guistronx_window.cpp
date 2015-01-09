#include "guistronx_window.h"
#include "ui_guistronx_window.h"

guistronx_window::guistronx_window(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::guistronx_window)
{
    ui->setupUi(this);
}

guistronx_window::~guistronx_window()
{
    delete ui;
}

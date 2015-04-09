#include "graphicsview.h"

class graphicsViewPrivate
{
public:
    graphicsViewPrivate()
    {
    }
};

graphicsView::graphicsView(QWidget *parent)
    : QGLWidget(QGLFormat(), parent)
    , d(new graphicsViewPrivate)
{
}

graphicsView::~graphicsView()
{
    delete d;
}

void graphicsView::initializeGL()
{
}

void graphicsView::resizeGL(int width, int height)
{
}

void graphicsView::paintGL()
{
}

#ifndef GRAPHICSVIEW_H
#define GRAPHICSVIEW_H

#include <QGLWidget>

class graphicsViewPrivate;

class graphicsView : public QGLWidget
{
    Q_OBJECT
public:
    explicit graphicsView(QWidget *parent = 0);
    ~graphicsView();

signals:

public slots:

protected:
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();

private:
    graphicsViewPrivate *d;
};

#endif // GRAPHICSVIEW_H

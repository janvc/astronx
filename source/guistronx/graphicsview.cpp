/*
 * Copyright 2015 Jan von Cosel
 *
 * This file is part of astronx.
 *
 * astronx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * astronx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have recieved a copy of the GNU General Public License
 * along with molconv. If not, see <http://www.gnu.org/licenses/>.
 *
 */


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

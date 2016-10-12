/************************************************************************
*
* Copyright 2014-2016 Michael McCrea (mtm5@uw.edu)
* Modification of QcScopeShm.h by Jakob Leben
*
* This file is part of Qt GUI for SuperCollider.
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*

* mtm - Modified to include code for FoaSpaceScope... work in progress.
* NOTE: for rollout, must also update:
  QtCollider/factories.cpp
  QtCollider/CMakeLists.txt
************************************************************************/

#ifndef QC_SCOPE_H
#define QC_SCOPE_H

#include "../QcHelper.h"

#include <SC_Types.h>

#include <QWidget>

// FIXME: Due to Qt bug #22829, moc can not process headers that include
// boost/type_traits/detail/has_binary_operator.hp from boost 1.48, so
// we have to wrap the shm interface into a separate class.

namespace QtCollider
{
  class ScopeShm;
}

using QtCollider::ScopeShm;

class QcSpaceScopeShm : public QWidget, QcHelper
{
  Q_OBJECT
  Q_PROPERTY( int serverPort READ serverPort WRITE setServerPort );
  Q_PROPERTY( int bufferNumber READ dummyInt WRITE setBufferNumber );
  Q_PROPERTY( float xOffset READ dummyFloat WRITE setXOffset );
  Q_PROPERTY( float yOffset READ dummyFloat WRITE setYOffset );
  Q_PROPERTY( float xZoom READ dummyFloat WRITE setXZoom );
  Q_PROPERTY( float yZoom READ dummyFloat WRITE setYZoom );
  Q_PROPERTY( int style READ style WRITE setStyle );
  Q_PROPERTY( QVariantList waveColors READ dummyVariantList WRITE setWaveColors );
  Q_PROPERTY( QColor background READ background WRITE setBackground );
  Q_PROPERTY( bool fill READ fill WRITE setFill );
  Q_PROPERTY( int updateInterval READ updateInterval WRITE setUpdateInterval );
  Q_PROPERTY( bool running READ running );
    // mtm
    Q_PROPERTY( int drawFrames READ drawFrames WRITE setDrawFrames );
    Q_PROPERTY( float ampShape READ ampShape WRITE setAmpShape );
    Q_PROPERTY( float normScale READ normScale WRITE setNormScale );
    Q_PROPERTY( float elWarp READ elWarp WRITE setElWarp );
    Q_PROPERTY( bool overlay READ overlay WRITE setOverlay );
    
    // TEMP
    Q_PROPERTY( float sat1 READ sat1 WRITE sat1);
    Q_PROPERTY( float val1 READ val1 WRITE val1);
    Q_PROPERTY( float al1 READ al1 WRITE al1);
    Q_PROPERTY( float al1max READ al1max WRITE al1max);
    
    Q_PROPERTY( float sat2 READ sat2 WRITE sat2);
    Q_PROPERTY( float val2 READ val2 WRITE val2);
    Q_PROPERTY( float al2 READ al2 WRITE al2);
    Q_PROPERTY( float al2max READ al2max WRITE al2max);
    Q_PROPERTY( float pos2 READ pos2 WRITE pos2);
    
    Q_PROPERTY( float sat3 READ sat3 WRITE sat3);
    Q_PROPERTY( float val3 READ val3 WRITE val3);
    Q_PROPERTY( float al3 READ al3 WRITE al3);
    Q_PROPERTY( float al3max READ al3max WRITE al3max);
    Q_PROPERTY( float pos3 READ pos3 WRITE pos3);

  public:
    QcSpaceScopeShm();
    ~QcSpaceScopeShm();
    int serverPort() const          { return _srvPort; }
    void setServerPort( int );
    void setBufferNumber( int );
    void setXOffset( float f )      { xOffset = f; }
    void setYOffset( float f )      { yOffset = f; }
    void setXZoom( float f )        { xZoom = f; }
    void setYZoom( float f )        { yZoom = f; }
//    int style() const               { return _style; }
//    void setStyle( int i )          { _style = i; }
    void setWaveColors( const QVariantList & colors );
    QColor background() const       { return _bkg; }
    void setBackground( const QColor &c ) { _bkg = c; update(); }
    bool fill()                     { return _fill; };
    void setFill(bool b)            { _fill = b; };
    int updateInterval() const;
    void setUpdateInterval( int i );
    bool running() const            { return _running; }
    QSize sizeHint() const          { return QSize( 500, 300 ); }
    QSize minimumSizeHint() const   { return QSize( 50, 50 ); }

    // mtm
    void setDrawFrames( int i )     { _drawFrames = i; }
    int drawFrames() const          { return _drawFrames; }
    void setAmpShape( float f )     { _ampShape = f; }
    float ampShape() const          { return _ampShape; }
    void setNormScale( float f )    { _normScale = f; }
    float normScale() const         { return _normScale; }
    void setElWarp( float f )       { _elWarp = f; updateElWarp(); }
    float elWarp() const            { return _elWarp; }
    void setOverlay( bool b )       { _overlay = b; resetPixmaps(); }
    bool overlay()                  { return _overlay; }
    int style() const               { return _style; }
    void setStyle( int i )          { _style = i; setMeterStyle(); }

    // TEMP
    void sat1( float f )     { _sat1 = f; }
    void val1( float f )     { _val1 = f; }
    void al1( float f )     { _al1 = f; }
    void al1max( float f )     { _al1max = f; }
    void sat2( float f )     { _sat2 = f; }
    void val2( float f )     { _val2 = f; }
    void al2( float f )     { _al2 = f; }
    void al2max( float f )     { _al2max = f; }
    void pos2( float f )     { _pos2 = f; }
    void sat3( float f )     { _sat3 = f; }
    void val3( float f )     { _val3 = f; }
    void al3( float f )     { _al3 = f; }
    void al3max( float f )     { _al3max = f; }
    void pos3( float f )     { _pos3 = f; }

    float sat1() const     { return _sat1; }
    float val1() const     { return _val1; }
    float al1() const     { return _al1; }
    float al1max() const     { return _al1max; }
    float sat2() const     { return _sat2; }
    float val2() const     { return _val2; }
    float al2() const     { return _al2; }
    float al2max() const     { return _al2max; }
    float pos2() const     { return _pos2; }
    float sat3() const     { return _sat3; }
    float val3() const     { return _val3; }
    float al3() const     { return _al3; }
    float al3max() const     { return _al3max; }
    float pos3() const     { return _pos3; }


  public Q_SLOTS:
    void start();
    void stop();
  protected:
    void resizeEvent ( QResizeEvent * );
    void paintEvent ( QPaintEvent * );

  private Q_SLOTS:
    void updateScope();
  private:
    void connectSharedMemory( int port );
    void initScopeReader( ScopeShm *, int index );
    void paint1D( bool overlapped, int channels, int maxFrames, int frames, const QRect &area, QPainter & );
    void paint2D( int channels, int maxFrames, int frames, const QRect &area, QPainter & );

    // mtm
    void lensDeform(const QPainterPath &source, const qreal radius, const qreal elev, const qreal directivity, const qreal lensDeform, QPainter &, bool fill);

    void paintAeda(int maxFrames, int availFrames, int frames, float ampShape, float normScale, bool overlay, bool fill, const QRect &area, QPainter &,
                   float sat1,
                   float val1,
                   float al1,
                   float al1max,
                   float sat2,
                   float val2,
                   float al2,
                   float al2max,
                   float pos2,
                   float sat3,
                   float val3,
                   float al3,
                   float al3max,
                   float pos3

                   );
    void paintAedaGrid( const QRect &area, QPainter & painter, bool overlay);

    void paintMercator(int maxFrames, int availFrames, int frames, float ampShape, float normScale, const QRect &area, QPainter & );
    void paintMercatorGrid(const QRect &area, QPainter & painter);

    void paintMollweide(int maxFrames, int availFrames, int frames, float ampShape, float normScale, const QRect &area, QPainter & );
    void paintMollweideGrid(const QRect &area, QPainter & painter);

    void prepPainter (QPainter &p, QColor &color, bool fill);
    void paintVectorScope(int maxFrames, int availFrames, int drawCount, const QRect &area, QPainter &painter);
    void paintVectorGrid(const QRect &area, QPainter & );

    // elevation warping for Mercator projection
    float mapCurvePi2 (float inVal, float curve);
    float mapCosPi2 (float inVal);
    float mapMercator(float inVal);

    void initMercatorPixDecoder( const int width, const int height, const float elWarp);
    void initMollweidePixDecoder( const int width, const int height);
    void storeCoeffs(int x, int y, float pat);
    void updateElWarp();
    void setMeterStyle();
    void resetPixmaps();

    int _drawFrames;
    float _ampShape, _normScale, _elWarp;
    bool _overlay;

    int nCols, nRows;
    bool decoderInitialized, gridDrawn;
    float wcoeff;
    float * xcoeffs, * ycoeffs, * zcoeffs;
    int storeDex, storeCnt;
    QImage img;
    QVector<QPoint> mwPixelPnts;
    
    // TEMP
    float  _sat1,
     _val1,
     _al1,
     _al1max,
     _sat2,
     _val2,
     _al2,
     _al2max,
     _pos2,
     _sat3,
     _val3,
     _al3,
     _al3max,
    _pos3;

    // end mtm

    int _srvPort;
    int _scopeIndex;

    ScopeShm *_shm;

    bool _running;
    float *_data;
    uint _availableFrames;

    QTimer *timer;
    float xOffset;
    float yOffset;
    float xZoom;
    float yZoom;
    int _style;
    QList<QColor> colors;
    QColor _bkg;
    bool _fill;

    QSize scopeSize;
    //mtm
//    QPixmap _pixmap;
    QPixmap _pixmap_meter;
    QPixmap _pixmap_bkg;
    QPixmap _pixmap_grid;
    QPixmap _prev_pixmap_meter;
};

#endif

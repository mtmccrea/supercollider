/************************************************************************
*
* Copyright 2014-2016 Michael McCrea (mtm5@uw.edu)
* Modification of QcScopeShm.cpp by Jakob Leben
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

#include "QcSpaceScopeShm.h"
#include "scope_shm_interface.hpp"
#include "../QcWidgetFactory.h"
#include "../debug.h"

#include <QPainter>
#include <QTimer>
#include <QResizeEvent>
#include <QWindow>
// mtm
//#include <QGenericMatrix>
#include <QGraphicsBlurEffect>
#include <QGraphicsScene>
#include <QGraphicsPixmapItem>

#define ROWPIX 120      /*  decoder width resolution     */
#define COLPIX 60       /*  decoder height resolution    */


QC_DECLARE_QWIDGET_FACTORY(QcSpaceScopeShm);

QcSpaceScopeShm::QcSpaceScopeShm() :
  _srvPort(-1),
  _scopeIndex(-1),
  _shm(new ScopeShm(this)),
  _running(false),
  _data(0),
  _availableFrames(0),
  xOffset( 0.f ),
  yOffset( 0.f ),
  xZoom( 1.f ),
  yZoom( 1.f ),
  _style( 0 ),
  _bkg( QColor(0,0,0) ),
//_fill( true ),

// mtm
  _fill( false ),
  _drawFrames( 4 ),
  _ampShape( 1.2 ),
  _normScale( 1.0),
  _elWarp( -1.0 ),
  _overlay( true ),

  nCols(1), nRows(1),
  decoderInitialized( false ),
  wcoeff(0.),
  xcoeffs(0), ycoeffs(0), zcoeffs(0),   // init to null pntr
  storeDex(0), storeCnt(0),
  img(),
  gridDrawn(false)

{
  setAttribute( Qt::WA_OpaquePaintEvent, true );
  setAutoFillBackground(false);

  setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );

  timer = new QTimer( this );
  timer->setInterval( 50 );
  connect( timer, SIGNAL( timeout() ), this, SLOT( updateScope() ) );
}

QcSpaceScopeShm::~QcSpaceScopeShm()
{
  stop();
}

void QcSpaceScopeShm::setServerPort( int port )
{
  if( _running ) {
    qcWarningMsg( "QScope: Can not change server port while running!" );
    return;
  }

  _srvPort = port;
}

void QcSpaceScopeShm::setBufferNumber( int n )
{
  if( _running ) {
    // TODO: release used reader?
    initScopeReader( _shm, n );
  }
  _scopeIndex = n;
}

void QcSpaceScopeShm::setWaveColors( const QVariantList & newColors )
{
  colors.clear();
  Q_FOREACH( const QVariant & var, newColors ) {
    QColor color = var.value<QColor>();
    if( !color.isValid() )
      colors.append( QColor( 0,0,0 ) );
    else
      colors.append( color );
  }
}

int QcSpaceScopeShm::updateInterval() const {
  return timer->interval();
}

void QcSpaceScopeShm::setUpdateInterval(int interval) {
  timer->setInterval( qMax(0, interval) );
}

void QcSpaceScopeShm::start()
{
    if( _running ) return;
    if( _srvPort < 0 || _scopeIndex < 0 ) return;

    connectSharedMemory( _srvPort );
    if( !_shm->client ) {
        stop();
        return;
    }

    initScopeReader( _shm, _scopeIndex );

    setMeterStyle(); // mtm

    timer->start();

    _running = true;
}

// this is called when setting meter style to initialize pixel decoders if their style is chosen
void QcSpaceScopeShm::setMeterStyle()
{
    switch (_style) {
        case 1:
            initMercatorPixDecoder( COLPIX * 2, COLPIX, _elWarp );
            break;
        case 2:
            initMollweidePixDecoder(COLPIX * 2, COLPIX);
            break;
    }

}

void QcSpaceScopeShm::stop()
{
  // TODO: release used reader?

  delete _shm->client;
  _shm->client = 0;

  timer->stop();

  _running = false;
}

void QcSpaceScopeShm::updateScope()
{
  bool valid = _shm->reader.valid();
  //qcDebugMsg(1, tr("valid = %1").arg(valid));
  if(!valid) return;

  bool ok = _shm->reader.pull( _availableFrames );
//   qcDebugMsg(1, tr("Got %1 frames").arg(_availableFrames) );
  if(ok) {
    _data = _shm->reader.data();
    update();
  }
}

void QcSpaceScopeShm::updateElWarp()
{
    initMercatorPixDecoder(COLPIX * 2, COLPIX, _elWarp);
}

void QcSpaceScopeShm::resetPixmaps()
{
    // TODO: do I need to cleanup former pixmaps? : _pixmap_bkg, _pixmap_grid, _prev_pixmap_meter
    _pixmap_meter = QPixmap(scopeSize);
    // create a black background
    _pixmap_bkg = QPixmap(scopeSize);
    _pixmap_bkg.fill(Qt::black);

    _prev_pixmap_meter = QPixmap(scopeSize);
    _prev_pixmap_meter.fill(Qt::transparent); // this happens in the paint event

    // create transparent pixmap for the grid/labels
    _pixmap_grid = QPixmap(scopeSize);
    _pixmap_grid.fill(Qt::transparent);

    gridDrawn = false; // trigger redrawing of the grid/label overlay
}

void QcSpaceScopeShm::resizeEvent ( QResizeEvent * ev )
{
    scopeSize = ev->size();

    // TODO: use pixmap.fill(_bkg) to handle trace, or remove
//    _pixmap = QPixmap(ev->size());
    resetPixmaps();
}


void QcSpaceScopeShm::paintEvent ( QPaintEvent * event )
{
    Q_UNUSED( event );

//    _pixmap.fill( _bkg ); // TODO: remove or replace with _pixmap_meter

    QPainter p;

    if( _running && _availableFrames ) {

        int chanCount = _shm->reader.channels();
        int maxFrames = _shm->reader.max_frames();

        QRect area (_pixmap_meter.rect());

        QPainter p_grid;
        QPainter p_meter;

        _prev_pixmap_meter = _pixmap_meter.copy();
        _pixmap_meter.fill(Qt::transparent);

        p_meter.begin(&_pixmap_meter);
        p_meter.setOpacity(_bkg.alphaF()); // trace amount set by user
        p_meter.drawPixmap(0,0,_prev_pixmap_meter);
        p_meter.setOpacity(1.0);

        // TODO: there may be a way to have this replace _prev_pixmap_meter
//        _pixmap_meter.fill( _bkg );

        // draw the meter display
        p_meter.begin(&_pixmap_meter);
        switch (_style) {
            case 0:
//                std::cout<<"painting case case\n";
                paintAeda(maxFrames, _availableFrames, _drawFrames, _ampShape, _normScale, _overlay, _fill, area, p_meter); break;
            case 1:
                paintMercator(maxFrames, _availableFrames, _drawFrames, _ampShape, _normScale, area, p_meter); break;
            case 2:
                paintMollweide(maxFrames, _availableFrames, _drawFrames, _ampShape, _normScale, area, p_meter); break;
            case 3:
                paintVectorScope(maxFrames, _availableFrames, _drawFrames, area, p_meter); break;
            case 4:
                paint1D(false, chanCount, maxFrames, _drawFrames, area, p_meter); break;

        }
        p_meter.end();

        // draw the overlay grid/labels
        if (!gridDrawn) {
            p_grid.begin(&_pixmap_grid);
            switch (_style) {
                case 0:
                    paintAedaGrid(area, p_grid, _overlay); break;
                case 1:
                    paintMercatorGrid(area, p_grid); break;
                case 2:
                    paintMollweideGrid(area, p_grid); break;
                case 3:
                    paintVectorGrid(area, p_grid); break;
            }
            p_grid.end();
            gridDrawn = true;
        }

        // paint the pixmap layers in order
        p.begin(this);
        p.drawPixmap(0, 0, _pixmap_bkg);
        p.drawPixmap(0, 0, _pixmap_meter);
        if(_style != 4) // no grid overlay on waveform
            p.drawPixmap(0, 0, _pixmap_grid);
    }
}

void QcSpaceScopeShm::paint1D( bool overlapped, int chanCount, int maxFrames, int frameCount,
                          const QRect &area, QPainter & painter )
{
    //qcDebugMsg( 0, tr("Drawing: data %1 / channels %2 / max-size %3").arg(_data!=0).arg(chanCount).arg(maxFrames) );
    if( frameCount < 2 || area.width() < 1 || area.height() < 1 ) return;

//std::cout<<"painting waveform\n";

    float yRatio = - yZoom * area.height() * 0.5;
    if( !overlapped ) yRatio /= chanCount;
    float yHeight = area.height();
    if( !overlapped ) yHeight /= chanCount;
    QPen pen;
    pen.setWidth(0);  // width==0 means width 1 regardless of transformations
    pen.setCapStyle( Qt::FlatCap );

    if( frameCount < area.width() )
    {
        float xRatio = xZoom * area.width() / (frameCount-1);

        for( int ch = 0; ch < chanCount; ch++ ) {
            float *frameData = _data + (ch * maxFrames); //frame vector
            float yOrigin = yHeight * (overlapped ? 0.5 : ch + 0.5);
            QColor strokeColor = ch < colors.count() ? colors[ch] : QColor(255,255,255);
            QColor fillColor(strokeColor); fillColor.setAlpha(0.65 * 255);
            pen.setColor( strokeColor );

            painter.save();
            painter.translate( area.x(), area.y() + yOrigin );
            painter.scale( xRatio, yRatio );
            painter.setPen(pen);

            QPainterPath path;
            path.moveTo( xOffset, frameData[0] );
            for( int f = 1; f < frameCount; ++f )
                path.lineTo( xOffset + f, frameData[f] );

            if (_fill) {
                path.lineTo(xOffset + frameCount, 0);
                path.lineTo(0, 0);
                path.lineTo(xOffset, frameData[0]);

                painter.fillPath(path, QBrush(fillColor));
            }

            painter.drawPath(path);

            painter.restore();
        }
    }
    else
    {
        int w = area.width();
        float fpp = frameCount / (float) w; // frames per x pixel
        float ypix = yRatio != 0.f ? -1/yRatio : 0.f; // value per y pixel;

        for( int ch = 0; ch < chanCount; ch++ )
        {
            float *frameData = _data + (ch * maxFrames); //frame vector
            float yOrigin = yHeight * (overlapped ? 0.5 : ch + 0.5);
            QColor strokeColor = ch < colors.count() ? colors[ch] : QColor(255,255,255);
            QColor fillColor(strokeColor); fillColor.setAlpha(0.65 * 255);

            painter.save();
            painter.translate( area.x(), area.y() + yOrigin );
            painter.setPen(pen);

            QPainterPath pathLine;
            QPainterPath pathFill;

            int p=0, f=1; // pixel, frame
            float min, max;
            min = max = frameData[0];

            while( p++ < w )
            {
                int f_max = fpp * p;

                for(; f < f_max; ++f)
                {
                    float d = frameData[f];
                    if( d < min ) min = d;
                    if( d > max ) max = d;
                }

                qreal x = p-1;
                float y = max * yRatio;
                pathLine.moveTo( x, y );
                y = qMax( min * yRatio, y+1 );
                pathLine.lineTo( x, y );

                if (_fill) {
                    pathFill.moveTo( x, y );
                    pathFill.lineTo( x, 0 );
                }

                // flip min/max to ensure continuity
                float val = min;
                min = max;
                max = val;
            }

            pen.setColor(strokeColor);
            painter.strokePath(pathLine, pen);

            if (_fill) {
                pen.setColor(fillColor);
                painter.strokePath(pathFill, pen);
            }

            painter.restore();
        }
    }
}

void QcSpaceScopeShm::paint2D( int chanCount, int maxFrames, int frameCount,
                          const QRect &area, QPainter & painter )
{
  int minSize = qMin( area.width(), area.height() );
  // NOTE: use yZoom for both axis, since both represent value, as opposed to index
  float xRatio = yZoom * minSize * 0.5;
  float yRatio = -yZoom * minSize * 0.5;
  QPoint center = area.center();
  QPen pen;
  pen.setWidth(0);  // width==0 means width 1 regardless of transformations
  pen.setCapStyle( Qt::FlatCap );
  pen.setColor(colors.count() ? colors[0] : QColor(255,255,255));

  painter.setPen(pen);
  painter.translate( center.x(), center.y() );
  painter.scale( xRatio, yRatio );

  QPainterPath path;

  if( chanCount >= 2 )
  {
    float *data1 = _data;
    float *data2 = _data + maxFrames;

    path.moveTo( data1[0], data2[0] );
    for( int f = 1; f < frameCount; ++f )
      path.lineTo( data1[f], data2[f] );
  }
  else
  {
    float *data1 = _data;
    path.moveTo( data1[0], 0.f );
    for( int f = 1; f < frameCount; ++f )
      path.lineTo( data1[f], 0.f );
  }

  painter.drawPath(path);
}

//mtm

typedef struct{ float mag,az,el,nrg; } maen;


void QcSpaceScopeShm::paintAedaGrid( const QRect &area, QPainter & painter, bool overlay )
{
    //std::cout<<"painting AEDA grid\n";

    float minSize, halfMin;
    if (overlay) {
        minSize = fmin(area.width(), area.height());
    } else {
        minSize = fmin(area.width()*0.5, area.height());
    }
    halfMin = minSize * 0.5;

    // we know we're in canvas bounds, so clipping=false is more efficient
    painter.setClipping(false);

    // center all drawing to top hemisphere view
    QPoint center = area.center();

    //TODO: could move a lot of this to an init funciton
    QString hemLabels[2] = {"Above","Below"};   // draws below first
    float el[] = { 0., 0.5236, 1.0472 };        // divide elevation pi/2 into thirds
    QString elDegs[] = {"0", "30", "60", "0", "-30", "-60"}; // draws below first

    // pixel spacing between hemispheres
    int   hemSpacing =          10;
    float hemSpacing_2 =        hemSpacing*0.5;
    float minSizeLessSpacing =  minSize - hemSpacing;
    float halfMinOffset =       -halfMin+hemSpacing_2;

    QRect txtRect(halfMinOffset, halfMinOffset, minSizeLessSpacing, minSizeLessSpacing);

    int gridcount;
    if (overlay) {
        gridcount = 1;
    } else {
        gridcount = 2;
    }

    for (int i=0; i<gridcount; i++) {
        float movey, movex;
        movey = center.y();

        if (overlay) {
            movex = center.x();
        } else {
            movex = center.x() - halfMin + (minSize*i); // draw right first
        }
        painter.translate(movex,  movey);


        float diam, rad, offset;
        QRectF bounds;

        //  Elevation Lines
        painter.setBrush(Qt::NoBrush);
        QColor gridColor(35,35,35); // very dark gray
        painter.setPen(gridColor);

        painter.rotate(-45);

        for ( int j=0; j<3; ++j) {
            diam   = (minSize-hemSpacing) * cos(el[j]);
            offset = -diam/2;
            bounds = QRectF(offset, offset, diam, diam);
            painter.drawEllipse(bounds);
            painter.drawText(bounds, Qt::AlignHCenter | Qt::AlignTop , elDegs[(i*3)+j]);
        };

        painter.rotate(45);

        /* draw azimuth angle labels */
        painter.setPen(QColor(75,75,75)); // light gray
        painter.drawText(txtRect, Qt::AlignHCenter | Qt::AlignTop ,   QString("0"));
        painter.drawText(txtRect, Qt::AlignHCenter | Qt::AlignBottom, QString("180"));
        painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignRight,  QString("-90"));
        painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignLeft,   QString("90"));

        // painter.setPen(Qt::red);
        if (!overlay)
            painter.drawText(txtRect, Qt::AlignLeft | Qt::AlignTop, hemLabels[i]);

        diam   = minSize-hemSpacing;
        rad    = diam/2;
        offset = -rad;
        bounds = QRectF(offset, offset, diam, diam);

//        QPoint fPnt( rad/5, rad/5 );   // focal point of gradient
        QPoint fPnt( 0, 0 );   // focal point of gradient
        if (i==1) fPnt *= -1;           // mirror across origin for lower hemi
        qreal colScale = 1 - 0.5*i;
        int col = colScale*255;
        QRadialGradient gr( 0, 0, rad, fPnt.x(), fPnt.y() );
        gr.setColorAt(0.0, QColor(col, col, col, 110));
        gr.setColorAt(0.4*colScale, QColor(col, col, col, 80));
        gr.setColorAt(0.9, QColor(col, col, col, 25));
        gr.setColorAt(1, QColor(0, 0, 0, 0));

        painter.setRenderHint(QPainter::Antialiasing);
        painter.setBrush(gr);
        painter.setPen(Qt::NoPen);
        painter.drawEllipse(bounds);

        // return to top left
        painter.translate(-movex,  -movey);
    }
}

void QcSpaceScopeShm::paintAeda( int maxFrames, int availFrames, int drawCount,
                           float ampShape, float normScale, bool overlay, bool fill,
                           const QRect &area, QPainter &painter )

{
    float minSize, halfMin;

    if (overlay) {
        minSize = fmin( area.width(), area.height() );
    } else {
        minSize = fmin( area.width() * 0.5, area.height() );
    }

    halfMin = minSize * 0.5;

    // we know we're in canvas bounds, so clipping=false is more efficient
//    painter.setClipping(false);

    // center all drawing to top hemisphere view
    QPoint center = area.center();

    // pixel spacing between hemispheres
    int   hemSpacing   = 10;
    float hemSpacing_2 = hemSpacing * 0.5;
//    float minSizeLessSpacing =  minSize - hemSpacing;
//    float halfMinOffset =       -halfMin+hemSpacing_2;
    float halfMinLessSpacing =  halfMin - hemSpacing_2;


    if (overlay) {
        painter.translate(center.x(), center.y());
    } else {
        // center on the left hemisphere
        painter.translate(center.x()-halfMin, center.y());
    }

    // painter.scale( halfMinLessSpacing, -halfMinLessSpacing  ); // scale coords 0 to 1

    // prepare for drawing aed circles
    painter.setPen(Qt::NoPen);
    QColor col;

    float *data1 = _data;
    float *data2 = _data + maxFrames;
    float *data3 = _data + (maxFrames * 2);
    float *data4 = _data + (maxFrames * 3);

    // clamp drawCount 1> availFrames
    int dcount = qMax(1, qMin(drawCount, availFrames));

    /* and even sampling of dcount samples withing the availFrames buffer */
//    //  an aeda array holding vectors of <a,e,d,amp>
//    aeda maenArray[dcount];
//
//    float frameHop = (float)availFrames/dcount;
//
//    /* get only dcount frames, evenly sampled within the data buffer */
//    // TODO: work out a max/min scheme instead, either max/min amp or directivity
//    for(int f = 0; f<dcount; ++f)
//    {
//        int dataDex = (int)(f * frameHop);
//        // data is now already aeda from Aeda ugen
//        maenArray[f].az  = data1[dataDex];
//        maenArray[f].el  = data2[dataDex];
//        maenArray[f].dir = data3[dataDex];
//        maenArray[f].amp = data4[dataDex];
//    }
//
//    /* sort in order of amplitude */
//    qsort(maenArray, dcount, sizeof(aeda), cmpfunc);
//
//    /* prepare normalize amps */
//    float max, sampAmp, normFac, normDist;
//    max = 0.0;
//    for (int i=0; i<dcount; ++i) {
//        sampAmp = maenArray[i].amp;
//        if (sampAmp > max) max = sampAmp;
//    }
//    normFac = 1. / max;
//    normDist = 1. - max;


    /* store max/min values of directivity */
    // TODO: sort by another parameter?

    // samples per resample window
    // force draw count to be even, rounded down
    int numMinMaxPairs = qMax( (int)(dcount * 0.5), 1);
    float winSamps = availFrames / numMinMaxPairs;
    dcount = numMinMaxPairs * 2;
    //  an aeda array holding vectors of <a,e,d,amp>
    maen maenArray[dcount];

    int drawDex, sampDex, minDex, maxDex, dex1, dex2;
    drawDex = sampDex = minDex = maxDex = dex1 = dex2 = 0;
    float min, max;

    for(int i=0; i<numMinMaxPairs; ++i) {
        int winStopDex;
        // initialize min and max frame
        // data1 is rE
        max = min = data1[sampDex];
        minDex = maxDex = sampDex;
        sampDex++;

        winStopDex = (int)(winSamps * (i+1));
//        minMaxWinSize = winStopDex-sampDex; // number of samples to iterate through for min/max

        for( ; sampDex<winStopDex; sampDex++) {
            float dir = data1[sampDex];
            if(dir<min) {min=dir; minDex=sampDex;};
            if(dir>max) {max=dir; maxDex=sampDex;};
//            sampDex++;
        };

//        for(int j=0; j<minMaxWinSize; ++j) {
//            float dir = data1[sampDex];
//            if(dir<min) {min=dir; minDex=sampDex;};
//            if(dir>max) {max=dir; maxDex=sampDex;};
//            sampDex++;
//        };

        if (maxDex<minDex) {
            dex1 = maxDex;
            dex2 = minDex;
        } else {
            dex1 = maxDex;
            dex2 = minDex;
        };

        maenArray[drawDex].mag = data1[dex1];
        maenArray[drawDex].az  = data2[dex1];
        maenArray[drawDex].el  = data3[dex1];
        maenArray[drawDex].nrg = data4[dex1];
        drawDex++;

        // could avoid confining to even dcount with:
//         if (drawDex == dcount) break;

        maenArray[drawDex].mag = data1[dex2];
        maenArray[drawDex].az  = data2[dex2];
        maenArray[drawDex].el  = data3[dex2];
        maenArray[drawDex].nrg = data4[dex2];
        drawDex++;
    }

// TODO: consider whether the aeda struct is necessary, now that not using qsort

    float xoff = minSize;
    float mag, a, e, nrg;

    for (int f = 0; f < dcount; ++f ) {
        /*  read in sorted aed data */
        mag    = maenArray[f].mag;
        a      = maenArray[f].az;
        e      = maenArray[f].el;
        nrg    = maenArray[f].nrg;

        // amplitude threshold, don't draw below -90 dB
        if(nrg > 0.00003) {

            /* save painter state, rotate coords to azimuth */
            painter.save();

            if (!overlay)
                if (e < 0.)
                    painter.translate(xoff,0);

            qreal diam   = minSize-hemSpacing;
            qreal rad    = diam/2;
            qreal offset = -rad;

            QPoint fPnt( 0, cos(e) * offset ); // focal point of gradient
            QRadialGradient gr(0, 0, rad, fPnt.x(), fPnt.y() );

            qreal mag_scale_inv = 1 - (mag*0.98);
            qreal mag_scale = mag * 0.99 + 0.01;
            qreal dcount_recip = 1 / sqrt(drawCount);

            //        std::cout << "mag_norm: " <<mag_norm<< "\tmag_scale: "<<mag_scale<< "\tmag_inv: "<<mag_scale_inv<<"\n";

            QColor col;
            int hue = 90 - (int)(30.*sin(e)); // 60 - yellow (positive elev), 120 - green (negative elev)

            col.setHsv(hue,
                       200*mag_scale_inv+55,   // saturation white > partial color
                       255,                 // value: "grayness"; black = 0
                       (90*mag_scale+55) *dcount_recip);   // alpha
            gr.setColorAt(0.0, col);

            col.setHsv(hue,
                       100*mag_scale_inv+155,    // saturation white > partial color
                       75*mag_scale_inv+180,    // value: "grayness"; black = 0
                       (40*mag_scale_inv+35) *dcount_recip);   // alpha
            gr.setColorAt(0.8 * pow(mag_scale_inv, 2) + 0.2, col);

            col.setHsv(hue,
                       255,       // saturation
                       125*mag_scale_inv+130, // grayness
                       (40*mag_scale_inv+25) *dcount_recip);  // alpha
            gr.setColorAt(0.6 * pow(mag_scale_inv, 2) + 0.35, col);

            col.setHsv(hue,
                       125,             // saturation, doesn't matter, it's black
                       0,               // grayness: black
                       0);   // alpha
//            gr.setColorAt(0.65 * mag_scale_inv + 0.3, col);
            gr.setColorAt(1.0, col);

//            gr.setColorAt(1.0, QColor(0, 0, 0, 0)); // black



//            qreal mag_scale_inv = 1 - (mag * 0.85);
//            qreal mag_scale = pow(mag, 2) * 0.9 + 0.1;
//            //        std::cout << "mag_norm: " <<mag_norm<< "\tmag_scale: "<<mag_scale<< "\tmag_inv: "<<mag_scale_inv<<"\n";
//
//            QColor col;
//            int hue = 90 - (int)(30.*sin(e)); // 60 - yellow (positive elev), 120 - green (negative elev)
//
//            col.setHsv(hue,
//                       mag_scale_inv*125,    // saturation white > partial color
//                       100*mag_scale+155,    // value: "grayness"
//                       130*mag_scale + 100); // alpha
//            gr.setColorAt(0.0, col);
//
//            col.setHsv(hue,
//                       130*mag_scale+125,    // saturation
//                       150*mag_scale+105,    // grayness
//                       100*mag_scale+60);    // alpha
//            gr.setColorAt(0.5 * mag_scale_inv, col);
//
//            col.setHsv(hue,
//                       125,                  // saturation, doesn't matter, it's black
//                       0,                    // grayness: black
//                       60*mag_scale );       // alpha
//            gr.setColorAt(1.0 * mag_scale_inv, col);
//
//            gr.setColorAt(1.0, QColor(0, 0, 0, 0)); // black


            //        painter.setRenderHint(QPainter::Antialiasing);
            painter.setBrush(gr);
            painter.setPen(Qt::NoPen);
            //        painter.setOpacity(1.0);

            // rotate negative azim to Ambisonics coords
            int aDeg = (int)(a/M_PI * -180.);
            painter.rotate(aDeg);

            QRectF bounds(offset, offset, diam, diam);

            painter.drawEllipse(bounds);

            // return to top left
            //        painter.translate(-movex,  -movey);

            painter.restore();
        }
    };
}


void QcSpaceScopeShm::paintVectorScope(int maxFrames, int availFrames, int drawCount,
                                  const QRect &area, QPainter &painter)
{
    float minSize = fmin( area.width(), area.height() );
    float halfMin = minSize * 0.5;

    // we know we're in canvas bounds, so clipping=false is more efficient
    //    painter.setClipping(false);

    // center all drawing to top hemisphere view
    QPoint center = area.center();

    // pixel spacing between hemispheres
    int   hemSpacing   = 10;
    float hemSpacing_2 = hemSpacing * 0.5;
//    float minSizeLessSpacing =  minSize - hemSpacing;
//    float halfMinOffset =       -halfMin+hemSpacing_2;
    float halfMinLessSpacing =  halfMin - hemSpacing_2;

    painter.translate(center.x(), center.y());

    painter.setBrush(Qt::yellow);

    float *data1 = _data;
    float *data2 = _data + maxFrames;

    int dcount = qMax(1, qMin(drawCount, availFrames)); // clamp drawCount 1> availFrames
    float frameHop = (float)availFrames/dcount;

    for (int f = 0; f < dcount; ++f ) {

        int dataDex = (int)(f * frameHop);
        float coh = data1[dataDex];
        float nrg = data2[dataDex];

        /* save painter state, rotate coords to azimuth */
        painter.save();

        qreal diam   = 5;      // diameter of plot point in pixels
        QPointF cenPnt( nrg, -coh);       // center of plot point, y inverted for screen coords
        cenPnt *= halfMinLessSpacing;     // scaled to pixel coords
        //        cenPnt += QPoint(offset, offset); // offset to corner of outlining ellipse
        QSizeF pntSize = QSizeF( diam, diam );
        QRectF pntRect ( cenPnt, pntSize );
        painter.setPen(Qt::yellow);
        painter.drawPoint(cenPnt);

        painter.restore();
    };
}

void QcSpaceScopeShm::paintVectorGrid( const QRect &area, QPainter &painter)
{
    float minSize, halfMin;
    minSize = fmin(area.width(), area.height());
    halfMin = minSize * 0.5;

    // we know we're in canvas bounds, so clipping=false is more efficient
    painter.setClipping(false);

    // center all drawing to top hemisphere view
    QPoint center = area.center();

    // pixel spacing between hemispheres
    int   hemSpacing =          10;
    float hemSpacing_2 =        hemSpacing*0.5;
    float minSizeLessSpacing =  minSize - hemSpacing;
    float halfMinOffset =       -halfMin+hemSpacing_2;

    float movey, movex;
    movey = center.y();
    movex = center.x();

    painter.translate(movex,  movey);

    //  Elevation Lines
    painter.setBrush(Qt::NoBrush);
    painter.setPen(Qt::red);

    QRectF bounds;
    QPen pen;
    pen.setColor(QColor(75,75,75,75));
    pen.setWidth(1);
    painter.setPen(pen);

        int numDivs = 8;
    int numDivLines = numDivs+1;

    float hatchStep = minSizeLessSpacing / numDivs;
    float hatchLen = hatchStep * 0.25;
    float halfHLen = hatchLen * 0.5;

    for ( int i=0; i<numDivLines; ++i) {
        // x lines
        int x = (int)(halfMinOffset + (hatchStep*i));
        painter.drawLine( x, -halfMinOffset, x, halfMinOffset);
        // y lines
        int y = (int)(halfMinOffset + (hatchStep*i));
        painter.drawLine( halfMinOffset, y, -halfMinOffset, y);

        // hatches
        float xs, xe, ys, ye;
        if (i==0) {xs = x; xe = xs+halfHLen;}
        else if (i == numDivs) {xs = x-halfHLen; xe = xs+halfHLen;}
        else {xs = x-halfHLen; xe = xs+hatchLen;}

        for ( int j=0; j<numDivLines; ++j) {
            y = (int)(halfMinOffset + (hatchStep*j));
            if (j==0) {ys = y; ye = ys+halfHLen;}
            else if (j == numDivs) {ys = y-halfHLen; ye = ys+halfHLen;}
            else {ys = y-halfHLen; ye = ys+hatchLen;}
            painter.drawLine( xs, y, xe, y); // horiz line
            painter.drawLine( x, ys, x, ye); // vert  line
        }
    }

    /* draw azimuth angle labels */
    painter.setPen(QColor(120,120,120));
    QRect txtRect(halfMinOffset, halfMinOffset, minSizeLessSpacing, minSizeLessSpacing);
    painter.drawText(txtRect, Qt::AlignHCenter | Qt::AlignTop ,   QString("Coherent"));
    painter.drawText(txtRect, Qt::AlignHCenter | Qt::AlignBottom, QString("Incoherent"));
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignRight,  QString("Pressure"));
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignLeft,   QString("Velocity"));

    // return to top left
    painter.translate(-movex,  -movey); // TODO: remove?
}

/*
    For each pixel (or pixel cluster), get the corresponding azimuth/elevation
    accourding to Mercator projection. Then calculated a cardioid decoder for that
    point in space and store the b-format coefficients.

    NOTE: pixel decoder coefficients need to be calculated only on init and resize
 */

/*
    Map 0>1 --> 0>pi/2
    Based on SC's CurveWarp, with an implicit spec range of 0>pi/2
    inVal is -1>1, sign is restored after mapping
 */
float QcSpaceScopeShm::mapCurvePi2(float inVal, float curve) {
    float absIn, a, grow, mappedVal;

    //  if curve is around 0, do linear warp
    if ((curve > -0.001) && (curve < 0.001)) {
        return (inVal / M_PI_2);
    }

    absIn = fabs(inVal);
    grow  = exp(curve);
    a     = M_PI_2 / (1-grow);                      // spec.range / (1.0 - grow);
    mappedVal = a - (a * pow(grow, absIn));
    return copysignf(mappedVal, inVal);         // restore the sign of the inVal
}

//float QcSpaceScopeShm::mapCosPi2(float inVal) {
//    float absIn, cosmap, linmap;
//
//    absIn = fabs(inVal);
//    cosmap = 0.5 - (cos(M_PI * absIn) * 0.5);
//    linmap = cosmap / M_PI_2;
//    return copysignf(linmap, inVal);            // restore the sign of the inVal
//}
//
////  Mercator "map" projection, i.e. sphere to cylinder
//float QcSpaceScopeShm::mapMercator(float inVal) {
//    float ypos = inVal * M_PI;  // yposition, -pi>pi
//    return asin(tanh( ypos ));  // Mercator result returns -pi/2>pi/2
//}


void QcSpaceScopeShm::storeCoeffs(int width, int height, const float pattern) {
    float xpos, ypos, theta, theta2, projAz, projEl;
    float projAzims[width];

    // TODO: refactor - store x/ypos terms outside the loop
    for (int j=0; j<width; ++j) {       // ... fill in the column coeffs
        xpos = j/(float)width;        // xposition, 0>1
        xpos = xpos*2 - 1;            // -1>1
        xpos = -xpos;                   // invert x for drawing coords
        projAzims[j] = xpos * M_PI;     // projected azimuth, partial term
    }

    for (int i=0; i <height; ++i) {     // for each row...
        ypos = i/(float)height;       // yposition, 0>1, top down
        ypos = ypos*2 - 1;            // yposition, -1>1
        ypos = -ypos;                   // invert y for drawing coords

        theta  = asin(ypos);
        theta2 = theta*2;

        projEl = asin((theta2 + sin(theta2)) / M_PI);   // projected elevation

        for (int j=0; j<width; ++j) {           // ... fill in the column coeffs
            projAz = projAzims[j] / cos(theta); // projected azimuth, complete the term

            // TODO: optimize this check
            if (((projAz <= M_PI)  && (projAz >= -M_PI)) &&
               ((projEl <= M_PI_2) && (projEl >= -M_PI_2))) {
                float costheta =  cos(projAz);
                float sintheta =  sin(projAz);
                float cosphipat = cos(projEl) * pattern;    // precomputed values for decoder coeffs
                float sinphipat = sin(projEl) * pattern;

                mwPixelPnts[storeDex] = QPoint(j, i);
                // index in 1D array of npoints size, rows first
                xcoeffs[storeDex] = cosphipat * costheta;
                ycoeffs[storeDex] = cosphipat * sintheta;
                zcoeffs[storeDex] = sinphipat;

                storeCnt++;
                storeDex++;
            }
        }
    }
}


/*  mono decoder coefficient calculation, use pattern = 0.5 for cardioid decode */
/*
    (1.0 - pattern) * 2.sqrt,        // W
    pattern * theta.cos * phi.cos,   // X
    pattern * theta.sin * phi.cos,   // Y
    pattern * phi.sin                // Z
*/

void QcSpaceScopeShm::initMollweidePixDecoder(const int width, const int height) {
    decoderInitialized = false;
    // delete any previously allocated memory for coefficients
    delete[] xcoeffs; delete[] ycoeffs; delete[] zcoeffs;

    // TODO: need to cleanup the previous QImage? QImageCleanupFunction?

    nCols = width;
    nRows = height;
    int npoints = nCols * nRows;
    img = QImage(nCols, nRows, QImage::Format_ARGB32) ;
    img.fill(Qt::black);

    // coefficients for decoder "points" are ordered in columns, bottom up
    // TODO: could use less points, Mollweide only uses ~73% of the total pixels

    const float pattern = 0.5;                      // decoder pattern, cardioid
    wcoeff  = M_SQRT1_2;                            // (1.0 - pattern) * sqrt(2.);
    xcoeffs = new (std::nothrow) float[npoints];    // new returns a pointer to head of array of npoint size
    ycoeffs = new (std::nothrow) float[npoints];
    zcoeffs = new (std::nothrow) float[npoints];

    // TODO: add error check:
    // if (xcoeffs == nullptr) { std::cout << "Error: memory could not be allocated for xcoeffs"; }
    // if (ycoeffs == nullptr) { std::cout << "Error: memory could not be allocated for ycoeffs"; }
    // if (zcoeffs == nullptr) { std::cout << "Error: memory could not be allocated for zcoeffs"; }

    // pixel coordinates of pixels to draw to
    mwPixelPnts = QVector<QPoint> (npoints);

    storeCnt = 0;
    storeDex = 0;

    storeCoeffs(nCols, nRows, pattern);

    //std::cout << "store count: "<<storeCnt<<"\n";

    decoderInitialized = true;
}


void QcSpaceScopeShm::initMercatorPixDecoder(const int width, const int height, const float elWarp) {
    //std::cout<<"initializing Mercator projection and allocating memory for coeffs \n";

    decoderInitialized = false;
    // delete any previously allocated memory for coefficients
    delete[] xcoeffs; delete[] ycoeffs; delete[] zcoeffs;

    // TODO: need to cleanup the previous QImage? QImageCleanupFunction?

    nCols = width;
    nRows = height;
    storeCnt = nCols * nRows;

    img = QImage(nCols, nRows, QImage::Format_ARGB32) ;

    wcoeff  = M_SQRT1_2;                            // (1.0 - pattern) * sqrt(2.);
    xcoeffs = new (std::nothrow) float[storeCnt];   // new returns a pointer to head of array of storeCnt size
    ycoeffs = new (std::nothrow) float[storeCnt];
    zcoeffs = new (std::nothrow) float[storeCnt];

    // TODO: add error check:
    // if (xcoeffs == nullptr) { std::cout << "Error: memory could not be allocated for xcoeffs"; }
    // if (ycoeffs == nullptr) { std::cout << "Error: memory could not be allocated for ycoeffs"; }
    // if (zcoeffs == nullptr) { std::cout << "Error: memory could not be allocated for zcoeffs"; }

    const float pattern = 0.5;                      // decoder pattern, cardioid
    float theta, phi;                               // az, el (longitude, lattitude)
    float costheta, sintheta;                       // precomputed trig ops
    float cosphipat, sinphipat;


    // store theta position for a single row, which is the same for all rows,
    // as well as pre-computed cos and sin of theta
    float thetas[nCols][3];

    // precompute for scaling each xpos 0>1 -> -1>1
    float scaleCols = 2.0 / nCols;

    for (int i=0; i<nCols; ++i) {
        float xpos = scaleCols * i - 1 ;    // xposition, -1>1
        float theta = -xpos * M_PI;         // -pi>pi, invert x so display left is positive

        thetas[i][0] = theta;
        thetas[i][1] = cos(theta);
        thetas[i][2] = sin(theta);
    }

    for (int i=0; i <nRows; ++i) {      // for each row...

        /* elevation/"lattitude"/yposition */
        float ypos = i/(float)nRows;     // yposition, 0>1, bottom>up
        ypos = ypos*2 - 1;                // yposition, -1>1
        ypos *= -1;                         // invert y because it will be packed top down in coefficient array
        // ypos *= M_PI;                       // yposition, -pi>pi, though Mercator result returns -pi/2>pi/2
        // phi = asin(tanh( ypos ));           // Mercator "map" projection, i.e. sphere to cylinder

        phi = mapCurvePi2(ypos, elWarp);      // warp to a curve
        //  phi = mapMercator(ypos);        // warp to Mercator projection
        //  phi = mapCosPi2(ypos);          // cosine warp

        cosphipat = cos(phi) * pattern;     // precomputed values for decoder coeffs
        sinphipat = sin(phi) * pattern;

        float rowhead = nCols*i;

        //  TODO: try a warp to the azimuth, where more width shown at 0 and rear is "shrunk"

        for (int j=0; j<nCols; ++j) {  // ... fill in the column coeffs
            theta =     thetas[j][0];
            costheta =  thetas[j][1];
            sintheta =  thetas[j][2];

            // index in 1D array of storeCnt size, rows first
            int dex = rowhead + j;
            xcoeffs[dex] = cosphipat * costheta; // *(xcoeffs + dex)
            ycoeffs[dex] = cosphipat * sintheta;
            zcoeffs[dex] = sinphipat;
        }
    }
    decoderInitialized = true;
}


void QcSpaceScopeShm::paintMollweide(int maxFrames, int availFrames,
                                int drawCount, float ampShape, float normScale,
                                const QRect &area, QPainter & painter )
{
    int minSize  = qMin(area.width()/2, area.height());
    int minSize2 = minSize*2;

    float *dataw =  _data;
    float *datax =  _data + maxFrames;
    float *datay =  _data + maxFrames*2;
    float *dataz =  _data + maxFrames*3;
    float *datapk = _data + maxFrames*4; // peak amplitude of energy

    // TODO: NOTE: currently just using first frame,
    // so "sample rate" is window refresh rate
    // It's therefore sensible to set input rms window to at least 1/windowRefreshRate
    float w,x,y,z, peak;
    w = dataw[0];
    x = datax[0];
    y = datay[0];
    z = dataz[0];
    peak = datapk[0]; // global max, slewed across frame on scsynth side

    // the w coefficient doesn't change with position, so no need for matrix
    float wdecoeff = wcoeff * w;

    QColor color = QColor(Qt::red);

    float decx, decy, decz, sum;
    float max = 0.0;
    float summedSamps[storeCnt]; // TODO: allocate on init?

    // decoded samples across the sphere
    for (int i=0; i<storeCnt; ++i) {
        decx = xcoeffs[i] * x;
        decy = ycoeffs[i] * y;
        decz = zcoeffs[i] * z;
        sum = decx + decy + decz + wdecoeff;
        summedSamps[i] = sum;
//        if(sum > max) max = sum;
    }
    // TODO: correct this measure of max amplitude of a decoded point
    max = sqrt((w * w * 2) + x*x + y*y + z*z) * 0.5;  // local max to this frame
    
//    float normDist =    1. - max;
    float fullNormFac = 1. / max;
    
    float refDist = peak - max;
    float normFac = (peak + (refDist * normScale)) / peak; // this tracks peak ("ref")
    
    for (int i=0; i<storeCnt; ++i) {
        float normVal, shapedVal;
        normVal = summedSamps[i] * fullNormFac;
        shapedVal = pow(normVal, ampShape);     // shape the amplitude response to accentuate peaks
        shapedVal *= max;                       // scale back to 0 > max
        shapedVal *= normFac;                   // scale from maxRef up to 1 according to normScale, maxRef = 1 when normFac=1

        color.setAlpha(shapedVal * 255);        // sum changes the opacity
        QPoint pnt = mwPixelPnts[i];
        img.setPixel(pnt.x(), pnt.y(), color.rgba());
    }
    

//    float normFac =  1. / max;
//    float normDist = 1. - max;
//
//    for (int i=0; i<storeCnt; ++i) {
//        float normVal, shapedVal, renormFac;
//        normVal = summedSamps[i] * normFac;
//        shapedVal = pow(normVal, ampShape);       // shape the amplitude response to accentuate peaks
//        renormFac = normScale * normDist + max;
//        shapedVal *= renormFac;
//        color.setAlpha(shapedVal * 255);    // sum changes the opacity
//        QPoint pnt = mwPixelPnts[i];
//        img.setPixel(pnt.x(), pnt.y(), color.rgba());
//    }

    /*added from init mercator*/

    // meter has its own scene with Pixmap (image) added to it
    QGraphicsScene meterscene;
    QGraphicsPixmapItem meter;
    // TODO: see Qt::NoFormatConversion (), http://doc.qt.io/qt-5/qt.html#GlobalColor-enum
    meter.setPixmap(QPixmap::fromImage(img));

    // blur the image for smooth upsizing
    // TODO: make blur radius a variable?
    int blurRad = 1;
    QGraphicsBlurEffect *blur = new QGraphicsBlurEffect;
    blur->setBlurRadius(blurRad);
    meter.setGraphicsEffect(blur);
    meterscene.addItem(&meter);

    // offset coordinates used to keep the meter centered on resize
    QPoint offsetPnt(area.width() - minSize2, area.height() - minSize);
    offsetPnt /= 2;  // center the scene in the view
    QSize scaleSize(minSize2, minSize); // add 2 px so boarder isn't cut off
    QRect offsetRect(offsetPnt, scaleSize);

    
    meterscene.render(&painter, offsetRect);

    painter.setPen(Qt::yellow);
    painter.drawText(offsetRect, Qt::AlignBottom | Qt::AlignLeft,
//                     QString::number(floor(max * pow(10., 4) + .5) / pow(10., 4)) + "\n" +
//                     QString::number(floor(peak * pow(10., 4) + .5) / pow(10., 4))+ "\n" +
                     QString::number(max) + "\n" +
                     QString::number(peak)+ "\n" +
                     QString::number(max-peak)
                     );
    /*end added from init mercator*/
}

void QcSpaceScopeShm::paintMollweideGrid( const QRect &area, QPainter &painter ) {

    //std::cout<<"painting MOLLWEIDE grid\n";

    int minSize  = qMin(area.width()/2, area.height());
    int minSize2 = minSize*2;

    // offset coordinates used to keep the meter centered on resize
    QPoint offsetPnt(area.width() - minSize2, area.height() - minSize);
    offsetPnt /= 2;  // center the scene in the view
    QSize scaleSize(minSize2, minSize); // add 2 px so boarder isn't cut off
    QRect offsetRect(offsetPnt, scaleSize);

    /* mask map edge with black ellipse */
    // maskWidth is as wide as a pixel row scaled up to new size
    int maskWidth = minSize2 * 1.5 / img.width() * 2;
    //    QColor maskColor (Qt::black);
    //    maskColor.setAlpha(125);
    //    QPen pen = QPen( maskColor );
    QPen pen(Qt::black);
    pen.setWidth(maskWidth);
    //    meterscene.addEllipse(maskWidth/2,maskWidth/2, meterscene.width(), meterscene.height(), pen);

    /* Mask Scene */
    QGraphicsScene maskScene;
    maskScene.setSceneRect(0, 0, minSize2, minSize);
    maskScene.addEllipse(maskWidth/5, maskWidth/5,
                         minSize2+maskWidth/5,
                         minSize+maskWidth/5,
                         pen );

    maskScene.render(&painter, offsetRect);
    /* Grid Scene */
    QGraphicsScene gridScene;
    gridScene.setSceneRect(0, 0,
                           minSize2 + maskWidth/5 - maskWidth,    // maskWidth = (2*sceneoffset),
                           minSize  + maskWidth/5 - maskWidth );

    float gSw   = gridScene.width();
    float gSw_2 = gridScene.width()/2.;
    float gSh   = gridScene.height();
    float gSh_2 = gridScene.height()/2.;
    int   sceneoffset = maskWidth/5 + maskWidth/2;


    // TODO: only draw once and call repaint on map item for update
    QColor gridColor = QColor(Qt::gray);
    gridColor.setAlpha(0.25*255);
    pen.setColor(gridColor);
    pen.setWidth(1);

    /* longitude lines */
    gridScene.addLine(gSw_2, 0, gSw_2, gSh, pen);   // longitude 0
    int numLongLines = 5;       // includes 0 and 180, though these aren't drawn in loop, so must be > 2
    for (int i=1; i<(numLongLines-1); ++i) {    // exclude long 0 & 180
        //    TODO: refactor
        float w = gSw * (float)i/(numLongLines-1);
        gridScene.addEllipse(gSw_2 - w/2., 0, w, gSh, pen);
    }

    /* lattitude lines */
    QBrush brush = QBrush(Qt::gray);
    // equator
    gridScene.addLine(0, gSh_2, gSw, gSh_2, pen);

    //  lattitude lines
    int numLatLines = 2;    // number of lattitude lines per hemisphere
    for (int j=0; j<numLatLines; ++j) {
        float el, theta, x, x2, y, y2;
        QPoint fmPnt, toPnt;

        el = (float)(j+1)/(numLatLines+1) * M_PI_2;   // ypos -pi/2>pi/2
        theta = asin(2*el/M_PI);

        x = cos(theta)*0.5; // -0.5>0.5, x = az*cos(theta) / PI; az = PI
        x2 = -x;
        y = sin(theta)*0.5;
        y2 = -y;
        x += 0.5;           // 0>1
        x2 += 0.5;
        y += 0.5;
        y2 += 0.5;
        x *= gSw;           // to scene coords
        x2 *= gSw;
        y *= gSh;
        y2 *= gSh;
        // draw lines at +/- elevation
        gridScene.addLine(x2, y, x, y, pen);
        gridScene.addLine(x2, y2, x, y2, pen);
    }

    // outline / longitude 180, fully opaque
    gridColor.setAlpha(125);
    pen.setColor(gridColor);
    gridScene.addEllipse(0, 0, gSw-1, gSh-1, pen); // -1 on w/h takes care of minor rounding errors in scene dimensions

    QPoint sceneOrigin(sceneoffset, sceneoffset);
    sceneOrigin += offsetPnt;

    gridScene.render( &painter, QRect(sceneOrigin, QSize( gSw, gSh )));

    // TODO: make graphic scene items so don't draw them every time
    painter.setPen(QColor(120,120,120)); // light gray
    QRect txtRect(sceneOrigin, QSize(gSw, gSh));
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignLeft, QString("ck"));
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignHCenter, QString("Front"));
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignRight, QString("Ba"));
    painter.drawText(txtRect, Qt::AlignTop     | Qt::AlignHCenter, QString("Above"));
    painter.drawText(txtRect, Qt::AlignBottom  | Qt::AlignHCenter, QString("Below"));

    txtRect = QRect(sceneOrigin, QSize(gSw_2, gSh));
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignHCenter, QString("Left"));

    txtRect = QRect(sceneOrigin + QPoint(gSw_2, 0), QSize(gSw_2, gSh));
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignHCenter, QString("Right"));
}


void QcSpaceScopeShm::paintMercator(int maxFrames, int availFrames,
                               int drawCount, float ampShape, float normScale,
                               const QRect &area, QPainter & painter ) {
    int minSize  = qMin(area.width()/2, area.height());
    int minSize2 = minSize*2;
    //    // we know we're in canvas bounds, so clipping=false is more efficient
    //    painter.setClipping(false);

    float summedSamps[storeCnt];    // TODO: allocate this in init?
    float *dataw = _data;
    float *datax = _data + maxFrames;
    float *datay = _data + (maxFrames * 2);
    float *dataz = _data + (maxFrames * 3);

    // NOTE: currently just using first frame,
    // so "sample rate" is window refresh rate
    // It's therefore sensible to set input rms window to at least 1/windowRefreshRate
    float w,x,y,z;
    w = dataw[0];
    x = datax[0];
    y = datay[0];
    z = dataz[0];
    // the w  coefficient doesn't change with position, so no need for matrix
    float wdecoeff = wcoeff * w;

    // TODO: check into QPainter::drawTiledPixmap
    // Both drawPixmap() and drawImage() produce the same result, except that drawPixmap() is faster on-screen

    QColor color = QColor(Qt::red);

    float decx, decy, decz, sum, max;
    for (int i=0; i<storeCnt; ++i) {
        decx = xcoeffs[i] * x;
        decy = ycoeffs[i] * y;
        decz = zcoeffs[i] * z;
        sum = decx + decy + decz + wdecoeff;
        summedSamps[i] = sum;
        if(sum > max) max = sum;
    }

    float normFac  = 1. / max;
    float normDist = 1. - max;

    for (int i=0; i<storeCnt; ++i) {
            float normVal, shapedVal, renormFac;
            normVal   = summedSamps[i] * normFac;
            shapedVal = pow(normVal, ampShape);       // shape the amplitude response to accentuate peaks
            renormFac = normScale * normDist + max;
            shapedVal *= renormFac;
            summedSamps[i] = shapedVal;
    }

    for (int row=0; row<nRows; ++row) {
        int offset = row*nCols;
        for (int col=0; col<nCols; ++col) {
            color.setAlpha(summedSamps[offset+col] * 255);    // sum changes the opacity
            img.setPixel(col, row, color.rgba());
        }
    }

    // meter has its own scene with Pixmap (image) added to it
    QGraphicsScene meterscene;
    QGraphicsPixmapItem meter;
    meter.setPixmap(QPixmap::fromImage(img));

    // blur the image for smooth upsizing
    // TODO: make blur radius a variable?
    int blurRad = 1;
    QGraphicsBlurEffect *blur = new QGraphicsBlurEffect;
    blur->setBlurRadius(blurRad);
    meter.setGraphicsEffect(blur);

    // add the meter to the scene
    meterscene.addItem(&meter);

    QPoint offsetPnt((area.width() - minSize2), area.height() - minSize);
    offsetPnt /= 2;  // center the scene in the view
    QSize scaleSize(minSize2, minSize);
    QRect offsetRect(offsetPnt, scaleSize);

    meterscene.render(&painter, offsetRect);
}

void QcSpaceScopeShm::paintMercatorGrid( const QRect &area, QPainter & painter ) {
    int minSize  = qMin(area.width()/2, area.height());
    int minSize2 = minSize*2;

    //std::cout<<"painting MERCATOR grid\n";

    QPoint offsetPnt((area.width() - minSize2), area.height() - minSize);
    offsetPnt /= 2;  // center the scene in the view
    QSize scaleSize(minSize2, minSize);
    QRect offsetRect(offsetPnt, scaleSize);

    // Grid Scene
    QGraphicsScene gridScene;
    gridScene.setSceneRect(QRect(QPoint(0,0), scaleSize));

    // TODO: only draw once and call repaint on map item for update
    QColor gridColor = QColor(Qt::gray);
    gridColor.setAlpha(0.25*255);
    QPen pen(gridColor);
    pen.setWidth(1);

    int numLongLines =  9;          // should be odd to include 0 longitude
    int numLatLines  =  5;          // should be odd to include equator
    int height = gridScene.height();
    int width  = gridScene.width();
    //  width-1 takes care of last line overshooting scene bounds
    float stepScale = (float)(width-1)/(numLongLines-1);

    for (int i=0; i<numLongLines; ++i) {
        int x = i * stepScale;
        gridScene.addLine(x, 0, x, height, pen);
    }

    stepScale = (float)(height-1)/(numLatLines-1);
    for (int i=0; i<numLatLines; ++i) {
        int y = i * stepScale;
        gridScene.addLine(0, y, width, y, pen);
    }

    //    meterscene.render(&painter, offsetRect);
    gridScene.render(&painter, offsetRect);

    // TODO: make graphic scene items so don't draw them every time
    painter.setPen(QColor(120,120,120)); // light gray
    QRect txtRect = offsetRect;
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignLeft,     QString("ck"));
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignHCenter,  QString("Front"));
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignRight,    QString("Ba"));
    painter.drawText(txtRect, Qt::AlignTop     | Qt::AlignHCenter,  QString("Above"));
    painter.drawText(txtRect, Qt::AlignBottom  | Qt::AlignHCenter,  QString("Below"));

    txtRect = QRect(offsetPnt, QSize(minSize,minSize));
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignHCenter,  QString("Left"));

    txtRect = QRect(offsetPnt.x()+minSize, offsetPnt.y(), minSize, minSize);
    painter.drawText(txtRect, Qt::AlignVCenter | Qt::AlignHCenter,  QString("Right"));
}

// end mtm

void QcSpaceScopeShm::connectSharedMemory( int port )
{
  try {
      server_shared_memory_client * client = new server_shared_memory_client(port);
      _shm->client = client;
      qcDebugMsg(1,"Shared memory connected");
  } catch (std::exception & e) {
      _shm->client = 0;
      qcErrorMsg(QStringLiteral("Cannot connect to shared memory: %1").arg(e.what()) );
  }
}

void QcSpaceScopeShm::initScopeReader( ScopeShm *shm, int index )
{
  shm->reader = shm->client->get_scope_buffer_reader( index );
  qcDebugMsg(1,QStringLiteral("Initialized scope buffer reader for index %1.").arg(index));
}



///* this method succeeds at clipping at the radius, but doesn't display negative elevation */
//void QcSpaceScopeShm::lensDeform(const QPainterPath &source, const qreal m_radius, const QRectF &rect, const qreal deformAmt, QPainter &painter)
//{
//    QPainterPath path;
//    path.addPath( source );
//
//    for (int i=0; i<path.elementCount(); ++i) {
//        const QPainterPath::Element &e = path.elementAt(i);
//
//        qreal x = e.x; // + offset.x();
//        qreal y = e.y; // + offset.y();
//
////        qreal dx = x; // - center.x(); // center is 0,0
////        qreal dy = y; // - center.y(); // center is 0,0
//        qreal hypot  = sqrt(x * x + y * y);
//        qreal len   = m_radius - hypot;
//        qreal ratio = len / m_radius;
//
////        qreal warpedx = x + (x * ratio);
////        qreal warpedy = y + (y * ratio);
//
////        path.setElementPositionAt(i, warpedx, warpedy);
//        qreal warpedx, warpedy;
//        if (len > 0) {
//            qreal ratio_scld = deformAmt * ratio;
//            warpedx = x + (x * ratio_scld);
//            warpedy = y + (y * ratio_scld);
//            path.setElementPositionAt(i, warpedx, warpedy);
//        } else {
//            // outside radius clips to radius
//            qreal clipScale = 1.-qAbs(len/hypot);
//            warpedx = x * clipScale;
//            warpedy = y * clipScale;
//            path.setElementPositionAt(i, warpedx, warpedy);
////            path.setElementPositionAt(i, x, y);
//        }
//
//    }
//
//    painter.drawPath(path);
//    //    return path;
//}

//// sort
//int cmpfunc (const void * a, const void * b)
//{
//    if ((*(maen*)a).nrg < (*(maen*)b).nrg)
//        return -1;
//    if ((*(maen*)a).nrg > (*(maen*)b).nrg)
//        return 1;
//    return 0;
//}


//  simplified from http://doc.qt.io/qt-5/qtwidgets-painting-deform-example.html
//  The center is 0,0 instead of mouse point, so dx and dy are just the x and y magnitudes

void QcSpaceScopeShm::lensDeform(const QPainterPath &source, const qreal radius,
                            const qreal elev, const qreal directivity, const qreal deformAmt,
                            QPainter &painter, bool fill) {

    QPainterPath primaryPath;
    primaryPath.addPath(source);
    QPainterPath foldedPath;
    foldedPath.addPath(source);

    bool folded = false;
    bool posElev = elev > 0;

    qreal radius_2 = radius*radius;

    for (int i=0; i<source.elementCount(); ++i) {
        const QPainterPath::Element &e = source.elementAt(i);

        qreal x     = e.x;
        qreal y     = e.y;
        qreal x_2   = x*x;
        qreal hypot = sqrt(x_2 + y*y);
        qreal len   = radius - hypot;

        qreal deformx, deformy;
        qreal ratio, ratio_scld, compressRatio;

        compressRatio = 1 - qAbs(len/hypot);


        if (len > 0) {                      /* source is within the radius */
            ratio = len / radius;
            ratio_scld = deformAmt*0.7 * ratio;
            deformx = x + (x * ratio_scld);
            deformy = y + (y * ratio_scld);
            // draw with lens deform
            primaryPath.setElementPositionAt(i, deformx, deformy);
            // expand to radius
            foldedPath.setElementPositionAt(i, x, -sqrt(radius_2 - x_2));
        } else {                            /* source is outside the radius */
            // projecting x/y out from radius, then folding back toward origin
            qreal yToRadius = sqrt(radius_2 - x_2);
            qreal yovershoot = qAbs(y) - yToRadius;
            qreal newHypot = radius - yovershoot;
            //            qreal comprRatio = newHypot / hypot;
            qreal comprRatio = newHypot / (radius + yovershoot);
            qreal newx = comprRatio * x;
            qreal newy = comprRatio * y;
            ratio = yovershoot / radius;

            //            //folding projection toward the origin
            //            qreal overshoot = hypot - radius;
            //            qreal newHypot = hypot - (2*overshoot);
            //            qreal scaleDim = newHypot / hypot;
            //            qreal newx = scaleDim * x;
            //            qreal newy = scaleDim * y;
            //
            ////            // folding projection parallel to y axis
            ////            qreal yToRadius = -sqrt(radius_2 - x_2);
            //////            qreal newy = y - (2*(y - yToRadius));
            ////            qreal newy = y + (2*(yToRadius - y));
            ////            qreal newHypot = sqrt(x_2 + newy*newy);
            ////            qreal scaleDim = radius / hypot;
            ////            qreal newx = scaleDim * x;
            //
            //            len = radius - newHypot;
            //            ratio = len / radius;

            ratio_scld = deformAmt * ratio;
            deformx = newx + (newx * ratio_scld);
            deformy = newy + (newy * ratio_scld);

            // draw with lens deform (folded in the boundary radius)
            foldedPath.setElementPositionAt(i, deformx, deformy);
            // compress to the radius boundary
            primaryPath.setElementPositionAt(i, x*compressRatio, y*compressRatio);

            folded = true;
        }
    }

    // set the color and painter's brush
    //    QColor posCol (170, 145, 57);
    //    QColor negCol (Qt::darkGreen);
    QColor posCol, negCol;
    int baseHue = 80;
    qreal dir = 230*pow(directivity, 7) + 25;
    qreal elScl = sin(elev)*45;
    posCol.setHsv( baseHue,
                  230,
                  200 + elScl,
                  dir );
    negCol.setHsv( (baseHue + (int)(0.167*255)) % 255,
                  230,
                  100 + elScl,
                  dir );

    // TODO: scale up to pixel coords here? So above calculations can be normalized

    // draw negative elevation first
    if (posElev) {
        // positive elevation is primary path
        if (fill) {
            // draw only positive response, filled
            prepPainter(painter, posCol, fill);
            painter.drawPath(primaryPath);
        } else {
            // not filled

            if (folded) {
                // there's a positive and negative response to draw
                prepPainter(painter, negCol, fill);     // negative will be drawn first
                painter.drawPath(foldedPath);           // so folded path is negative
                prepPainter(painter, posCol, fill);     // and so drawn first
                painter.drawPath(primaryPath);
            } else {
                // positive only
                prepPainter(painter, posCol, fill);     // and so drawn first
                painter.drawPath(primaryPath);
            }
        }
    } else {
        // negative elevation is primary path
        prepPainter(painter, negCol, fill);
        painter.drawPath(primaryPath);
        if (folded) {
            prepPainter(painter, posCol, fill); // and so drawn first
            painter.drawPath(foldedPath);
        }
    }
}

void QcSpaceScopeShm::prepPainter ( QPainter &p, QColor &color, bool fill) {
    if (fill) {
        p.setPen(Qt::NoPen);
        p.setBrush(color);
    } else {
        p.setBrush( Qt::NoBrush );
        p.setPen(color);
    }
}

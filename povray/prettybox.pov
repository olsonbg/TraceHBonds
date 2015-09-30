#version 3.6;
//
// Edit the bottom of this file so that it includes the appropriate
// file generated from TraceHBonds program.
//
// povray +W640 +H480 +A0.1 +AM2 -D prettybox.pov +ONPT200K-640x480.png
// povray +W640 +H480 -D prettybox-animate.ini +O200K/Frame-640x480-
//
// Add +UA for transparent backgound
//------------------------------------------------------------------------
#include "colors.inc"
#include "textures.inc"
#include "shapes.inc"
#include "transforms.inc"

//------------------------------------------------------------------------
#macro Camera_LookAt( Xlook, Ylook, Zlook )
camera {
      angle 14
      location  <Xlook/2,Ylook/2,-900>
      right     x*image_width/image_height
//    Negative 'right' makes it a right handed coordinate system.
//    with +z coming out of the screen.
      look_at   <Xlook/2,Ylook/2,Zlook/2>
      Rotate_Around_Trans(<20.0,clock*360,0>, <Xlook/2,Ylook/2,Zlook/2>)
}
#end

//------------------------------------------------------------------------
// sun -------------------------------------------------------------------
light_source{ <67.4504/2,67.4504/2,-900>
              color White
              Rotate_Around_Trans(<20.0,clock*360,0>,<67.4504/2,67.4504/2,67.4417/2>)
            }
// sky -------------------------------------------------------------------
sky_sphere{ pigment{ gradient <0,1,0>
                     color_map{ [0   color rgb<1,1,1>         ]//White
                                [0.4 color rgb<0.14,0.14,0.46>]//~Navy
                                [0.6 color rgb<0.14,0.14,0.46>]//~Navy
                                [1.0 color rgb<1,1,1>         ]//White
                              }
                     scale 3 }
           } // end of sky_sphere


#declare ChainRadius = 0.7;

#declare ChainInvisiblePigment = pigment{color rgbt<0,0,0,1>};

#declare ChainRGB =
 array[10]
 {
  rgb<0.2317689248,0.2983710996,0.7530599677>,
  rgb<0.3605146877,0.4890730716,0.9057955249>,
  rgb<0.5093102463,0.6574314477,0.9917652973>,
  rgb<0.6660533453,0.78655854,1.0>,
  rgb<0.8119437768,0.8612552278,0.9340117644>,
  rgb<0.9271841987,0.8706813455,0.8003321552>,
  rgb<0.9941061637,0.8097466668,0.6169175421>,
  rgb<0.9990261307,0.679438159,0.4039675231>,
  rgb<0.93355859,0.484634481,0.1756267707>,
  rgb<0.7957405835,0.21340006,0.0>,
 }

#declare ChainLength3pigment  = pigment{color ChainRGB[ 0]};
#declare ChainLength5pigment  = pigment{color ChainRGB[ 1]};
#declare ChainLength7pigment  = pigment{color ChainRGB[ 2]};
#declare ChainLength9pigment  = pigment{color ChainRGB[ 3]};
#declare ChainLength11pigment = pigment{color ChainRGB[ 4]};
#declare ChainLength13pigment = pigment{color ChainRGB[ 5]};
#declare ChainLength15pigment = pigment{color ChainRGB[ 6]};
#declare ChainLength17pigment = pigment{color ChainRGB[ 7]};
#declare ChainLength19pigment = pigment{color ChainRGB[ 8]};
#declare ChainLength21pigment = pigment{color ChainRGB[ 9]};
#declare ChainLength23pigment = pigment{color ChainRGB[ 9]};

#declare ChainFinish  = finish {ambient 0.15 diffuse 0.85 phong 1};

// #declare ChainLength3  = texture{ pigment{ChainInvisiblePigment}};
#declare ChainLength3  = texture{ finish {ChainFinish} pigment{ChainLength3pigment}  };
#declare ChainLength5  = texture{ finish {ChainFinish} pigment{ChainLength5pigment}  };
#declare ChainLength7  = texture{ finish {ChainFinish} pigment{ChainLength7pigment}  };
#declare ChainLength9  = texture{ finish {ChainFinish} pigment{ChainLength9pigment}  };
#declare ChainLength11 = texture{ finish {ChainFinish} pigment{ChainLength11pigment} };
#declare ChainLength13 = texture{ finish {ChainFinish} pigment{ChainLength13pigment} };
#declare ChainLength15 = texture{ finish {ChainFinish} pigment{ChainLength15pigment} };
#declare ChainLength17 = texture{ finish {ChainFinish} pigment{ChainLength17pigment} };
#declare ChainLength19 = texture{ finish {ChainFinish} pigment{ChainLength19pigment} };
#declare ChainLength21 = texture{ finish {ChainFinish} pigment{ChainLength21pigment} };
#declare ChainLength23 = texture{ finish {ChainFinish} pigment{ChainLength23pigment} };

#macro PBC( Xlen, Ylen, Zlen )
object {
	// Wire_Box(A, B, WireRadius, Merge)
	Wire_Box(<0,0,0>,<Xlen, Ylen, Zlen>, 0.2, 0)
	texture{
		pigment{ color rgb<0.75,0.55,>}
		finish {  phong 1}
	}
	rotate<0,0,0>
	translate<0,0,0>
}
#end // end PBC macro

#macro Legend( Xcoord, Ycoord, Zcoord )
	#declare Yoffset=5;
	#declare Xoffset=3;

	#declare TextureArray = array[11] {ChainLength3,  ChainLength5,
	                                   ChainLength7,  ChainLength9,
	                                   ChainLength11, ChainLength13,
	                                   ChainLength15, ChainLength17,
	                                   ChainLength19, ChainLength21,
	                                   ChainLength23 }

	#declare Index = 0;
	#while ( Index < 11 )
		sphere_sweep {
			linear_spline
			3,
			<  Xcoord          , Ycoord, Zcoord>, 0.7
			<  Xcoord+Xoffset  , Ycoord, Zcoord>, 0.7
			<  Xcoord+2*Xoffset, Ycoord, Zcoord>, 0.7
			tolerance 0.07
			texture{TextureArray[Index]}
		}

		#declare Ycoord=Ycoord-Yoffset;
		#declare Index=Index+1;
	#end // of while
#end // Legend macro

//Legend( 60.0, 50.0, 0.0)
#include "HBonds1.pov"

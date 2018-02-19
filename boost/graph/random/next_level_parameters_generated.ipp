// Copyright (C) 2011-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine

#include <math.h>
void step_coeffs_internal(double a,double b,double c,double d,double db,
double crea_par[4])
{
  double t1;
  double t103;
  double t104;
  double t107;
  double t108;
  double t11;
  double t111;
  double t112;
  double t113;
  double t117;
  double t120;
  double t123;
  double t124;
  double t13;
  double t132;
  double t133;
  double t134;
  double t135;
  double t136;
  double t138;
  double t139;
  double t140;
  double t141;
  double t142;
  double t143;
  double t144;
  double t145;
  double t146;
  double t147;
  double t148;
  double t149;
  double t15;
  double t150;
  double t152;
  double t153;
  double t154;
  double t155;
  double t156;
  double t157;
  double t158;
  double t159;
  double t16;
  double t160;
  double t162;
  double t163;
  double t166;
  double t167;
  double t168;
  double t169;
  double t170;
  double t171;
  double t172;
  double t174;
  double t175;
  double t176;
  double t177;
  double t178;
  double t18;
  double t180;
  double t181;
  double t182;
  double t184;
  double t185;
  double t188;
  double t19;
  double t190;
  double t191;
  double t192;
  double t194;
  double t195;
  double t199;
  double t2;
  double t201;
  double t202;
  double t203;
  double t205;
  double t209;
  double t21;
  double t210;
  double t212;
  double t213;
  double t215;
  double t217;
  double t219;
  double t22;
  double t220;
  double t222;
  double t226;
  double t231;
  double t232;
  double t233;
  double t235;
  double t236;
  double t237;
  double t238;
  double t239;
  double t24;
  double t240;
  double t245;
  double t247;
  double t248;
  double t250;
  double t252;
  double t253;
  double t254;
  double t255;
  double t256;
  double t257;
  double t258;
  double t259;
  double t26;
  double t260;
  double t261;
  double t262;
  double t263;
  double t264;
  double t265;
  double t27;
  double t270;
  double t271;
  double t272;
  double t273;
  double t274;
  double t277;
  double t279;
  double t280;
  double t281;
  double t286;
  double t287;
  double t288;
  double t29;
  double t291;
  double t294;
  double t296;
  double t297;
  double t298;
  double t30;
  double t300;
  double t303;
  double t306;
  double t307;
  double t308;
  double t31;
  double t310;
  double t311;
  double t312;
  double t314;
  double t317;
  double t318;
  double t32;
  double t320;
  double t322;
  double t326;
  double t33;
  double t333;
  double t337;
  double t338;
  double t339;
  double t341;
  double t342;
  double t344;
  double t348;
  double t350;
  double t352;
  double t36;
  double t360;
  double t363;
  double t365;
  double t366;
  double t367;
  double t369;
  double t37;
  double t370;
  double t371;
  double t376;
  double t378;
  double t379;
  double t380;
  double t381;
  double t386;
  double t390;
  double t392;
  double t394;
  double t4;
  double t40;
  double t403;
  double t407;
  double t409;
  double t410;
  double t415;
  double t416;
  double t418;
  double t420;
  double t422;
  double t425;
  double t426;
  double t43;
  double t431;
  double t433;
  double t434;
  double t435;
  double t440;
  double t443;
  double t444;
  double t449;
  double t45;
  double t452;
  double t464;
  double t47;
  double t479;
  double t488;
  double t49;
  double t492;
  double t494;
  double t496;
  double t498;
  double t5;
  double t502;
  double t504;
  double t506;
  double t51;
  double t516;
  double t52;
  double t520;
  double t521;
  double t523;
  double t525;
  double t529;
  double t531;
  double t533;
  double t54;
  double t541;
  double t545;
  double t548;
  double t549;
  double t551;
  double t552;
  double t553;
  double t558;
  double t56;
  double t561;
  double t562;
  double t565;
  double t568;
  double t57;
  double t572;
  double t575;
  double t578;
  double t59;
  double t592;
  double t6;
  double t60;
  double t607;
  double t61;
  double t612;
  double t617;
  double t62;
  double t621;
  double t63;
  double t633;
  double t65;
  double t650;
  double t654;
  double t66;
  double t664;
  double t666;
  double t668;
  double t673;
  double t675;
  double t677;
  double t68;
  double t686;
  double t70;
  double t708;
  double t72;
  double t721;
  double t723;
  double t724;
  double t729;
  double t73;
  double t732;
  double t735;
  double t738;
  double t74;
  double t75;
  double t750;
  double t76;
  double t77;
  double t78;
  double t79;
  double t8;
  double t80;
  double t809;
  double t81;
  double t82;
  double t83;
  double t831;
  double t840;
  double t845;
  double t846;
  double t851;
  double t86;
  double t87;
  double t88;
  double t89;
  double t9;
  double t92;
  double t93;
  double t94;
  double t95;
  double t96;
  double t99;
  {
    t1 = a*db;
    t2 = d*db;
    t4 = log(c+a+t1+b+d+t2);
    t5 = t4/8.0;
    t6 = b*db;
    t8 = log(b+a+t1+d+t6+c+t2);
    t9 = t8/8.0;
    t11 = log(c+a+t1+b+d);
    t13 = c*db;
    t15 = log(b+a+t1+t13+t6+c+d);
    t16 = t15/8.0;
    t18 = log(c+a+t1+b+t13+d);
    t19 = t18/8.0;
    t21 = log(t2+a+t1+b+c+t13+d);
    t22 = t21/8.0;
    t24 = a+b+c+d;
    t26 = log((1.0+db)*t24);
    t27 = t26/8.0;
    t29 = log(b+a+t1+d+t6+c);
    t30 = t29/8.0;
    t31 = t26-t21-t15-t11+t4+t18-t8+t29;
    t32 = t31/2.0;
    t33 = 1/db;
    t36 = db*db;
    t37 = 1/t36;
    t40 = 1/t36/db;
    t43 = log(a+b+t2+c+t13+d);
    t45 = log(a+c+b+t6+t2+t13+d);
    t47 = log(a+b+c+d+t2);
    t49 = log(a+b+c+t13+d);
    t51 = log(a+t13+b+t6+c+d);
    t52 = log(t24);
    t54 = log(a+d+b+t6+c+t2);
    t56 = log(a+d+b+t6+c);
    t57 = t26+t43-t45-t47-t49+t51+t29+t52-t21+t18+t54-t56+t4-t11-t15-t8;
    t59 = t36*t36;
    t60 = 1/t59;
    t61 = t57*t60/8.0;
    t62 = t5-t9-t11/8.0-t16+t19-t22+t27+t30+t32*t33+3.0/4.0*t31*t37+t32*t40+t61
;
    t63 = 1/d;
    t65 = 1/c;
    t66 = 1/b;
    t68 = a*a;
    t70 = t65*t66*t68*a;
    t72 = t15/3.0;
    t73 = t26/3.0;
    t74 = t29/3.0;
    t75 = t8/3.0;
    t76 = 4.0/3.0*t26;
    t77 = t18/3.0;
    t78 = t4/3.0;
    t79 = t11/3.0;
    t80 = 4.0/3.0*t8;
    t81 = 4.0/3.0*t15;
    t82 = 4.0/3.0*t29;
    t83 = t21/3.0;
    t86 = 2.0*t15;
    t87 = 2.0*t26;
    t88 = 2.0*t8;
    t89 = 2.0*t29;
    t92 = t51/3.0;
    t93 = t56/3.0;
    t94 = t54/3.0;
    t95 = t45/3.0;
    t96 = t92-t11+t18-t21+t76-t93-t80-t81+t82+t94-t95+t4;
    t99 = t57*t60/3.0;
    t103 = 4.0/3.0*t18;
    t104 = 4.0/3.0*t21;
    t107 = 2.0*t21;
    t108 = 2.0*t18;
    t111 = t49/3.0;
    t112 = t43/3.0;
    t113 = t76-t111+t92-t104-t95-t8-t11+t29-t81+t112+t103+t4;
    t117 = 4.0/3.0*t4;
    t120 = 2.0*t4;
    t123 = t47/3.0;
    t124 = -t15-t104+t76-t123+t117+t18+t112+t94-t11+t29-t80-t95;
    t132 = t15/4.0;
    t133 = t8/4.0;
    t134 = t26/4.0;
    t135 = t29/4.0;
    t136 = t26+t29-t15-t8;
    t138 = t4/4.0;
    t139 = t21/4.0;
    t140 = t54/4.0;
    t141 = 3.0/2.0*t8;
    t142 = t45/4.0;
    t143 = 3.0/2.0*t15;
    t144 = t51/4.0;
    t145 = t18/4.0;
    t146 = t11/4.0;
    t147 = t56/4.0;
    t148 = 3.0/2.0*t26;
    t149 = 3.0/2.0*t29;
    t150 = t138-t139+t140-t141-t142-t143+t144+t145-t146-t147+t148+t149;
    t152 = t51/2.0;
    t153 = t21/2.0;
    t154 = t11/2.0;
    t155 = t56/2.0;
    t156 = t4/2.0;
    t157 = t18/2.0;
    t158 = t54/2.0;
    t159 = t45/2.0;
    t160 = t152+t26-t8-t153-t154-t155-t15+t156+t157+t158-t159+t29;
    t162 = t57/4.0;
    t163 = t162*t60;
    t166 = t65*b;
    t167 = (-t132-t133+t134+t135+t136*t33+t150*t37+t160*t40+t163)*t63*t166;
    t168 = t15/2.0;
    t169 = t26/2.0;
    t170 = t29/2.0;
    t171 = t8/2.0;
    t172 = t157+t87+t170-t171-t153-t86;
    t174 = 3.0*t15;
    t175 = 3.0/2.0*t21;
    t176 = 3.0*t26;
    t177 = 3.0/2.0*t18;
    t178 = -t174+t156-t159-t175+t176+t177-t154+t152-t141+t149;
    t180 = t49/2.0;
    t181 = t43/2.0;
    t182 = -t180-t11+t87-t175-t45+t177-t141-t86-t155+t51+t181+t4+t149+t158;
    t184 = t57/2.0;
    t185 = t184*t60;
    t188 = t156-t168-t88-t153+t170+t87;
    t190 = 3.0/2.0*t4;
    t191 = 3.0*t8;
    t192 = -t175-t159+t157+t176-t154+t149+t190-t143+t158-t191;
    t194 = t47/2.0;
    t195 = -t194-t155+t152+t149+t190-t175+t54-t143-t11-t88+t181-t45+t87+t18;
    t199 = -t15+t18-t21+t26;
    t201 = t43/4.0;
    t202 = t49/4.0;
    t203 = -t143+t177+t201-t175-t202-t133+t148+t138-t142-t146+t144+t135;
    t205 = -t180-t159+t26-t21-t154-t171+t18+t181+t152+t170-t15+t156;
    t209 = (t145-t132-t139+t134+t199*t33+t203*t37+t205*t40+t163)*t63*c;
    t210 = t157+t87-t171+t156-t168-t107;
    t212 = 3.0*t21;
    t213 = -t212+t176+t181-t159-t143+t190-t141+t170+t177-t154;
    t215 = -t180+t87+t152+t190-t11-t194-t143+t177+t43-t141-t45-t107+t158+t29;
    t217 = -t21+t26+t4-t8;
    t219 = t47/4.0;
    t220 = t145+t148+t201-t146-t219+t190-t141+t140-t132-t142-t175+t135;
    t222 = t4-t154-t8+t157-t194-t159+t26+t170-t168+t158+t181-t21;
    t226 = (-t139+t134+t138-t133+t217*t33+t220*t37+t222*t40+t163)*d*t65;
    t231 = t45/24.0;
    t232 = t15/24.0;
    t233 = t54/24.0;
    t235 = t51/24.0;
    t236 = t29/24.0;
    t237 = t26/24.0;
    t238 = t8/24.0;
    t239 = -t26+t15+t45+t8-t51-t54+t56-t29;
    t240 = t239/6.0;
    t245 = -t57;
    t247 = t245*t60/24.0;
    t248 = t231+t232-t233+t56/24.0-t235-t236-t237+t238+t240*t33+t239*t37/4.0+
t240*t40+t247;
    t250 = b*b;
    t252 = t65*t250*b;
    t253 = t248*t63*t252;
    t254 = t45/6.0;
    t255 = t51/6.0;
    t256 = t15/6.0;
    t257 = t26/6.0;
    t258 = 2.0/3.0*t26;
    t259 = t56/6.0;
    t260 = 2.0/3.0*t51;
    t261 = t54/6.0;
    t262 = t8/6.0;
    t263 = 2.0/3.0*t45;
    t264 = t29/6.0;
    t265 = 2.0/3.0*t15;
    t270 = t21/6.0;
    t271 = t49/6.0;
    t272 = t43/6.0;
    t273 = t18/6.0;
    t274 = t263+t270-t158+t171-t258+t271+t155-t260-t272-t273+t265-t170;
    t277 = t245*t60/6.0;
    t279 = (t254-t255+t256-t257+(-t258+t259-t260-t261+t262+t263-t264+t265)*t33+
(t155-t158-t26+t15+t171+t45-t170-t51)*t37+t274*t40+t277)*t63;
    t280 = 2.0/3.0*t54;
    t281 = 2.0/3.0*t8;
    t286 = t4/6.0;
    t287 = t47/6.0;
    t288 = -t152-t258+t168+t281+t263+t270-t286-t280-t272-t170+t287+t155;
    t291 = (t254+t262-t261-t257+(-t280-t258-t255+t256+t281+t259-t264+t263)*t33+
(t8-t170-t54+t45-t152-t26+t155+t168)*t37+t288*t40+t277)*t65;
    t294 = -t51-t26+t15+t45;
    t296 = 3.0/2.0*t45;
    t297 = 3.0/2.0*t51;
    t298 = t139-t140+t296-t148+t133+t143+t202-t135-t145-t297+t147-t201;
    t300 = t180-t158+t45-t26+t15+t153+t155-t157-t181+t171-t51-t170;
    t303 = -t162*t60;
    t306 = (t132-t144-t134+t142+t294*t33+t298*t37+t300*t40+t303)*t63*c;
    t307 = 2.0*t45;
    t308 = t171+t168-t152-t158-t87+t307;
    t310 = 3.0*t45;
    t311 = 3.0/2.0*t54;
    t312 = t153-t170-t176+t310+t155-t181-t311-t297+t141+t143;
    t314 = t307-t311-t87+t180+t141-t297+t143-t156-t157+t56-t29-t43+t21+t194;
    t317 = -t184*t60;
    t318 = t45+t8-t26-t54;
    t320 = t139-t144+t296+t141-t148+t147-t201+t132+t219-t138-t135-t311;
    t322 = t168-t26-t152+t8+t194+t45-t156-t181-t170+t155-t54+t153;
    t326 = (t133-t134-t140+t142+t318*t33+t320*t37+t322*t40+t303)*d*t65;
    t333 = t180+t263-t258-t261+t153-t260-t264-t181+t262+t259+t265-t157;
    t337 = c*c;
    t338 = (-t257-t255+t256+t254+(t270+t265-t260+t271-t272+t263-t258-t273)*t33+
(t45-t157-t181+t180+t153-t51-t26+t15)*t37+t333*t40+t277)*t63*t337;
    t339 = -t87+t153-t152+t168+t307-t181;
    t341 = 3.0/2.0*t43;
    t342 = -t176+t175-t341+t310-t158+t171+t143-t297-t157+t180;
    t344 = -t156+t143+t49+t307-t87-t341-t54-t18-t297+t175-t170+t155+t194+t8;
    t348 = -t158+t171+t153-t181-t87+t307;
    t350 = t310-t176+t175+t141-t341+t168+t194-t156-t152-t311;
    t352 = t15-t87-t157+t180+t307+t155-t51-t341+t141-t170-t4+t47-t311+t175;
    t360 = t256-t258-t255+t281-t156+t263-t264-t181+t153+t194+t259-t280;
    t363 = d*d;
    t365 = (t254+t262-t261-t257+(t281-t286-t258+t263+t270-t280+t287-t272)*t33+(
t153+t8-t181+t45-t26-t156+t194-t54)*t37+t360*t40+t277)*t363*t65;
    t366 = t21/24.0;
    t367 = t18/24.0;
    t369 = t43/24.0;
    t370 = t45+t49+t21+t15-t43-t51-t26-t18;
    t371 = t370/6.0;
    t376 = t231+t366-t367-t235+t232+t49/24.0-t237-t369+t371*t33+t370*t37/4.0+
t371*t40+t247;
    t378 = t337*c;
    t379 = t376*t63*t378;
    t380 = 2.0/3.0*t21;
    t381 = 2.0/3.0*t43;
    t386 = t180-t261-t286-t258+t168+t262+t263-t157+t380-t381-t152+t287;
    t390 = t21-t26-t43+t45;
    t392 = t133-t144+t175-t138+t202+t296-t145-t148-t140-t341+t132+t219;
    t394 = -t26+t45+t168+t171-t156-t152-t43-t157-t158+t180+t194+t21;
    t403 = t271+t263-t258+t171-t156-t158-t273+t256-t381+t194-t255+t380;
    t407 = t4/24.0;
    t409 = t8-t43-t54-t26+t47-t4+t45+t21;
    t410 = t409/6.0;
    t415 = t366+t238+t231-t407-t233+t47/24.0-t237-t369+t410*t33+t409*t37/4.0+
t410*t40+t247;
    t416 = t363*d;
    t418 = t415*t416*t65;
    t420 = (t379+(t270-t257-t272+t254+(t380+t271-t255-t258+t263-t381+t256-t273)
*t33+(-t43+t45+t21-t26-t152+t180+t168-t157)*t37+t386*t40+t277)*t337+(-t134+t139
-t201+t142+t390*t33+t392*t37+t394*t40+t303)*d*c+(t270-t257-t272+t254+(t263+t262
-t286-t381+t380-t261+t287-t258)*t33+(t171-t156+t45-t26-t158+t21+t194-t43)*t37+
t403*t40+t277)*t363+t418)*t66;
    t422 = 1/a;
    t425 = -t31;
    t426 = t425/6.0;
    t431 = t238+t366-t236+t11/24.0-t367+t232-t407-t237+t426*t33+t425*t37/4.0+
t426*t40+t247;
    t433 = t431*t63*t70;
    t434 = 2.0/3.0*t18;
    t435 = t11/6.0;
    t440 = -t258+t254-t255-t434+t380-t272+t171-t156+t154+t265-t170+t271;
    t443 = (t270-t273+t256-t257+(-t434+t435+t380-t258-t286+t262+t265-t264)*t33+
(t154-t26+t15+t171-t170+t21-t18-t156)*t37+t440*t40+t277)*t63;
    t444 = 2.0/3.0*t4;
    t449 = -t261-t170+t287-t258+t380-t157+t254-t272-t444+t154+t168+t281;
    t452 = (-t286-t257+t262+t270+(t281-t273+t435-t258-t444-t264+t380+t256)*t33+
(t154+t21-t170-t157+t8+t168-t26-t4)*t37+t449*t40+t277)*t65;
    t464 = (t139-t145-t134+t132-t199*t33-t203*t37-t205*t40+t303)*t63*c;
    t479 = (-t138+t139+t133-t134-t217*t33-t220*t37-t222*t40+t303)*d*t65;
    t488 = t82-t83+t77+t54-t45-t79+t51-t80-t81+t76+t78-t56;
    t492 = t65*t250;
    t494 = -t86-t159+t87+t152+t170-t171;
    t496 = t157-t296-t174+t149-t153-t141+t297-t155+t158+t176;
    t498 = t149+t18+t54-t180-t141-t21+t181+t156-t296+t87+t297-t86-t56-t154;
    t502 = -t159+t170+t87+t158-t88-t168;
    t504 = t149-t153+t156-t143+t311+t176-t191-t296-t155+t152;
    t506 = t157+t149-t56-t143-t194+t87+t181-t21-t296-t88+t311+t51+t4-t154;
    t516 = -t264-t286-t181+t380+t435-t258-t434+t262+t159+t180-t152+t265;
    t520 = (t256-t257+t270-t273+(-t258+t254-t434-t255+t380+t271+t265-t272)*t33+
(-t18+t159-t181-t26-t152+t21+t15+t180)*t37+t516*t40+t277)*t63*t337;
    t521 = -t87-t181-t157+t168+t159+t107;
    t523 = t171-t176-t156-t177-t341+t212-t152+t296+t180+t143;
    t525 = -t177-t170+t107+t296+t194-t51-t158+t49-t341-t87+t8-t4+t154+t143;
    t529 = -t87+t171+t107+t159-t181-t156;
    t531 = t296-t176+t194+t141-t157-t190+t168-t341+t212-t158;
    t533 = -t18-t152+t47-t87+t107-t170+t296-t190-t341+t154-t54+t141+t180+t15;
    t541 = t435-t273+t194+t380-t181+t256+t159+t281-t258-t444-t158-t264;
    t545 = (-t257-t286+t270+t262+(-t444-t261+t281+t380-t258+t287-t272+t254)*t33
+(-t4-t158-t26+t21+t194+t8-t181+t159)*t37+t541*t40+t277)*t363*t65;
    t548 = t51/8.0;
    t549 = t54/8.0;
    t551 = t45/8.0;
    t552 = -t239;
    t553 = t552/2.0;
    t558 = t548+t549+t30-t9-t56/8.0+t27-t551-t16+t553*t33+3.0/4.0*t552*t37+t553
*t40+t61;
    t561 = 4.0/3.0*t45;
    t562 = 4.0/3.0*t51;
    t565 = 2.0*t51;
    t568 = t112+t29+t77-t81+t76+t562+t54-t561-t56-t83-t111-t8;
    t572 = 4.0/3.0*t54;
    t575 = 2.0*t54;
    t578 = -t123+t76-t83+t29+t78+t51+t572-t561-t56-t80-t15+t112;
    t592 = (-t132-t142+t144+t134-t294*t33-t298*t37-t300*t40+t163)*t63*c;
    t607 = (-t133-t142+t134+t140-t318*t33-t320*t37-t322*t40+t163)*d*t65;
    t612 = 2.0/3.0*t29;
    t617 = t254+t265-t255+t281-t156+t154+t259-t612+t153-t261-t258-t157;
    t621 = (-t257+t262-t264+t256+(t281+t270-t273-t258-t612+t265+t435-t286)*t33+
(-t26+t154+t153+t8-t29-t157-t156+t15)*t37+t617*t40+t277)*t63*t65;
    t633 = (-t134+t132+t133-t135-t136*t33-t150*t37-t160*t40+t303)*t63*t166;
    t650 = t281-t273-t152+t435-t158-t258+t155-t286+t270+t265-t612+t159;
    t654 = (-t257+t262-t264+t256+(-t261+t254+t259-t258-t612+t265-t255+t281)*t33
+(t8+t159-t29+t155+t15-t26-t152-t158)*t37+t650*t40+t277)*t63*t492;
    t664 = -t153-t159+t87-t86+t152+t157;
    t666 = -t296-t180-t171+t177+t176+t297+t170-t175-t174+t181;
    t668 = -t8+t297-t49-t86-t155+t87-t154+t29+t43+t158+t177+t156-t296-t175;
    t673 = t88-t87+t153-t158-t156+t159;
    t675 = -t170+t194+t168-t181+t296-t311+t191-t190-t176+t175;
    t677 = t88-t190+t155-t43-t152+t296-t87-t29+t47-t311+t175+t154+t15-t157;
    t686 = t76+t51-t49-t81-t75-t79+t103+t74+t43-t104+t78-t45;
    t708 = -t49-t81+t562+t76-t75+t94-t93+t43-t21+t74+t18-t561;
    t721 = t43/8.0;
    t723 = -t370;
    t724 = t723/2.0;
    t729 = -t22+t721-t551+t27+t19+t548-t49/8.0-t16+t724*t33+3.0/4.0*t723*t37+
t724*t40+t61;
    t732 = 4.0/3.0*t43;
    t735 = 2.0*t43;
    t738 = t51-t75+t76-t49+t732-t15+t18-t104-t123+t94+t78-t561;
    t750 = (-t142-t139+t134+t201-t390*t33-t392*t37-t394*t40+t163)*d*c;
    t809 = -t104-t72+t76-t47+t74+t117-t45+t43-t79+t77+t54-t80;
    t831 = -t21-t47+t74+t76+t4-t93-t72+t43+t572+t92-t561-t80;
    t840 = t76+t77+t4+t92-t561-t8-t47+t732-t104-t72+t54-t111;
    t845 = -t409;
    t846 = t845/2.0;
    t851 = -t551+t5-t47/8.0+t549+t721-t22+t27-t9+t846*t33+3.0/4.0*t845*t37+t846
*t40+t61;
    crea_par[0] = 1.0/4.0+t62*t63*t70+((-t72+t73+t74-t75+(t76+t77+t78-t79-t80-
t81+t82-t83)*t33+(t4-t86+t87+t18-t88-t11-t21+t89)*t37+t96*t40+t99)*t63*t65+((-
t72+t77-t83+t73+(t74+t103-t79-t104-t75+t76-t81+t78)*t33+(-t86-t8-t107+t108-t11+
t87+t4+t29)*t37+t113*t40+t99)*t63+(t73-t75-t83+t78+(t117+t77-t79-t80-t72+t76-
t104+t74)*t33+(t120+t18+t29+t87-t107-t88-t15-t11)*t37+t124*t40+t99)*t65)*t66)*
t68+(t167+(-t168+t169+t172*t33+t178*t37+t182*t40+t185)*t63+(-t171+t169+t188*t33
+t192*t37+t195*t40+t185)*t65+(t209-t153+t169+t210*t33+t213*t37+t215*t40+t185+
t226)*t66)*a+(t253+(t279+t291)*t250+(t306-t169+t159+t308*t33+t312*t37+t314*t40+
t317+t326)*b+t338+(-t169+t159+t339*t33+t342*t37+t344*t40+t317)*c+(-t169+t159+
t348*t33+t350*t37+t352*t40+t317)*d+t365+t420)*t422;
    crea_par[1] = 1.0/4.0+t433+(t443+t452)*t66*t68+(t167+(t464-t169+t153-t210*
t33-t213*t37-t215*t40+t317+t479)*t66)*a+(t74-t75+t73-t72+(-t81-t93+t76-t95+t94+
t82+t92-t80)*t33+(-t45+t89-t56+t51-t88-t86+t87+t54)*t37+t488*t40+t99)*t63*t492+
((-t168+t169+t494*t33+t496*t37+t498*t40+t185)*t63+(-t171+t169+t502*t33+t504*t37
+t506*t40+t185)*t65)*b+(t520+(t153-t169+t521*t33+t523*t37+t525*t40+t317)*c+(
t153-t169+t529*t33+t531*t37+t533*t40+t317)*d+t545)*t66+(t558*t63*t252+((t92+t73
-t95-t72+(-t81+t94+t74-t561-t93+t76-t75+t562)*t33+(-t307+t565+t29-t86-t56-t8+
t87+t54)*t37+t568*t40+t99)*t63+(-t75+t73+t94-t95+(t92+t76-t561+t572-t93+t74-t72
-t80)*t33+(t29-t15-t307-t88+t51-t56+t87+t575)*t37+t578*t40+t99)*t65)*t250+(t592
+t169-t159-t308*t33-t312*t37-t314*t40+t185+t607)*b+t420)*t422;
    crea_par[2] = 1.0/4.0+t433+(t621+t452*t66)*t68+(t633+(t171-t169-t188*t33-
t192*t37-t195*t40+t317)*t65+(t209+t479)*t66)*a+t654+(t171-t169-t502*t33-t504*
t37-t506*t40+t317)*t65*b+(t169-t168+t664*t33+t666*t37+t668*t40+t185)*t63*c+(
t171-t169+t673*t33+t675*t37+t677*t40+t317)*d*t65+((t73+t77-t72-t83+(t112-t104+
t76-t95+t92+t103-t81-t111)*t33+(-t49-t107+t87+t108-t45+t51-t86+t43)*t37+t686*
t40+t99)*t63*t337+(-t153+t169-t521*t33-t523*t37-t525*t40+t185)*c+t545)*t66+(
t253+t291*t250+(t592+t326)*b+(-t95+t92-t72+t73+(-t111+t112+t77-t561+t76-t81-t83
+t562)*t33+(-t49-t21+t43+t87-t86+t565-t307+t18)*t37+t708*t40+t99)*t63*t337+(-
t159+t169-t339*t33-t342*t37-t344*t40+t185)*c+t365+(t729*t63*t378+(t73+t112-t95-
t83+(-t104-t561+t76-t111-t72+t92+t732+t77)*t33+(-t307-t107+t735-t15+t87+t51+t18
-t49)*t37+t738*t40+t99)*t337+t750+t418)*t66)*t422;
    crea_par[3] = 1.0/4.0+t433+(t621+t443*t66)*t68+(t633+(-t169+t168-t172*t33-
t178*t37-t182*t40+t317)*t63+(t464+t226)*t66)*a+t654+(-t169+t168-t494*t33-t496*
t37-t498*t40+t317)*t63*b+(-t169+t168-t664*t33-t666*t37-t668*t40+t317)*t63*c+(-
t171+t169-t673*t33-t675*t37-t677*t40+t185)*d*t65+(t520+(t169-t153-t529*t33-t531
*t37-t533*t40+t185)*d+(-t83+t73-t75+t78+(t76-t123+t112+t94-t80-t95-t104+t117)*
t33+(t54-t107+t120+t43-t45+t87-t47-t88)*t37+t809*t40+t99)*t363*t65)*t66+(t253+
t279*t250+(t306+t607)*b+t338+(t169-t159-t348*t33-t350*t37-t352*t40+t185)*d+(t73
+t94-t95-t75+(t76-t123+t78-t83-t80+t112-t561+t572)*t33+(t43-t47-t307+t87+t4-t88
+t575-t21)*t37+t831*t40+t99)*t363*t65+(t379+t750+(t73+t112-t83-t95+(-t104-t561+
t76-t123-t75+t94+t78+t732)*t33+(t54-t8+t735-t307+t4-t107+t87-t47)*t37+t840*t40+
t99)*t363+t851*t416*t65)*t66)*t422;
    return;
  }
}


#include <stdio.h>
#include <math.h>
#include <SDL.h>
#include "pga3d.h"
typedef struct vert{
    float x, y, z;
} vert;
typedef struct cube{
    vert v[8];
}cube;
int main()
{
    SDL_Window *w = NULL;
    SDL_Renderer *r = NULL;
    SDL_Event  e;
    uint32_t dt, t;
    cube c;
    int i;
    double a;
    c.v[0] = (vert) {300, 500, 500};
    c.v[1] = (vert) {500, 500, 500};
    c.v[2] = (vert) {500, 300, 500};
    c.v[3] = (vert) {300, 300, 500};
    
    c.v[4] = (vert) {300, 500, 300};
    c.v[5] = (vert) {500, 500, 300};
    c.v[6] = (vert) {500, 300, 300};
    c.v[7] = (vert) {300, 300, 300};
    
    w = SDL_CreateWindow(
        "spinny",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        800, 800,
        SDL_WINDOW_RESIZABLE
    );
    if (!w) goto exit;

    r = SDL_CreateRenderer(
        w, -1 , SDL_RENDERER_ACCELERATED
    );
    if (!r) goto exit;
    t = SDL_GetTicks();
    a = 0;
    do{
        SDL_SetRenderDrawColor(r, 0, 0, 0, 0xff);
        SDL_RenderClear(r);

        SDL_SetRenderDrawColor(r, 0, 0xFF, 0x80, 0xff);

        for (i =0; i < 4; i++){
            SDL_RenderDrawLine(r,
                c.v[i].x, c.v[i].y,
                c.v[(i+1)%4].x, c.v[(i+1)%4].y
            );
        }
        SDL_SetRenderDrawColor(r, 0xFF, 0, 0x80, 0xff);

        for (i = 0; i < 8; i++)
        {
            multivector_t mv = point(
                c.v[i].x,
                c.v[i].y,
                c.v[i].z
            );

            multivector_t rotor = PGA3D_rotor(0.01, PGA3D_wedge(e2(), e3()));
            rotor = PGA3D_mul(rotor, PGA3D_rotor(0.01, PGA3D_wedge(e3(), e1())));
            mv = PGA3D_mul(rotor, PGA3D_mul( mv, PGA3D_reverse(rotor)));

            c.v[i].x = mv.mvec[13];
            c.v[i].y = mv.mvec[12];
            c.v[i].z = mv.mvec[11];
        }
        
        for (i =0; i < 4; i++){
            SDL_RenderDrawLine(r,
                c.v[i+4].x, c.v[i+4].y,
                c.v[(i+1)%4+4].x, c.v[(i+1)%4+4].y
            );
        }       
        
        SDL_SetRenderDrawColor(r, 0xFF, 0xDD, 0x30, 0xff);
        for (i =0; i < 4; i++){
            SDL_RenderDrawLine(r,
                c.v[i].x, c.v[i].y,
                c.v[i+4].x, c.v[i+4].y
            );
        }       
        SDL_RenderPresent(r);

        SDL_PollEvent(&e);
        SDL_Delay(2);
        dt = SDL_GetTicks() - t;
        t = SDL_GetTicks();
    }
    while(e.type != SDL_QUIT);
    
exit:
    SDL_DestroyRenderer(r);
    SDL_DestroyWindow(w);
    return 0;
}

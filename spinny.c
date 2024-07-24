#include <stdio.h>
#include <string.h>
#include <math.h>
#include <SDL.h>
#include "pga3d.h"
multivector_t spring(multivector_t M, multivector_t B, multivector_t anchor, multivector_t corner, float k);
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
    SDL_Rect anchor_box;
    multivector_t M, B, anchor_mv, dM, dB;

    c.v[0] = (vert) {-0.25, 0.25, 0.25};
    c.v[1] = (vert) {0.25, 0.25, 0.25};
    c.v[2] = (vert) {0.25, -0.25, 0.25};
    c.v[3] = (vert) {-0.25, -0.25, 0.25};
    
    c.v[4] = (vert) {-0.25, 0.25, -0.25};
    c.v[5] = (vert) {0.25, 0.25, -0.25};
    c.v[6] = (vert) {0.25, -0.25, -0.25};
    c.v[7] = (vert) {-0.25, -0.25, -0.25};
    
    anchor_box.w = 8;
    anchor_box.h = 8;
    anchor_box.x = 400-4;
    anchor_box.y = 400-4;

    anchor_mv = PGA3D_point( 0, 0.25, 0);

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
    dt = 0;
    memset(&M, 0, sizeof(M));
    M.mvec[0] =1;
    B = PGA3D_add(
        PGA3D_wedge(PGA3D_e1(), PGA3D_e2()),
        PGA3D_smul( 1.3,
                PGA3D_wedge(PGA3D_e1(), PGA3D_e3())
        )
    );
    

    do{
        SDL_SetRenderDrawColor(r, 0, 0, 0, 0xff);
        SDL_RenderClear(r);

        SDL_SetRenderDrawColor(r, 0, 0xFF, 0x10, 0xff);

        SDL_RenderDrawRect(r, &anchor_box);
        SDL_SetRenderDrawColor(r, 0, 0xFF, 0x80, 0xff);

        for (i =0; i < 4; i++){
            SDL_RenderDrawLine(r,
                (c.v[i].x +1)* 400, (c.v[i].y +1)* 400,
                (c.v[(i+1)%4].x +1)* 400, (c.v[(i+1)%4].y +1)* 400
            );
        }
        SDL_SetRenderDrawColor(r, 0xFF, 0, 0x80, 0xff);

        dM = PGA3D_smul(
            -0.9,
            PGA3D_mul(M, B)
        );

        dB = PGA3D_sub(
            spring(M, B,
                anchor_mv,
                PGA3D_point(
                    c.v[0].x,
                    c.v[0].y,
                    c.v[0].z
                ),
            12
            ), 
            PGA3D_smul(0.5,
                PGA3D_dual( 
                    PGA3D_reverse(
                    PGA3D_sub(
                        PGA3D_mul(PGA3D_dual(B), B),
                        PGA3D_mul(B, PGA3D_dual(B))
                    )
                )
                )
            )
        );

        M = PGA3D_add(M, PGA3D_smul(dt/1000.0, dM));
        B = PGA3D_add(B, PGA3D_smul(dt/1000.0, dB));
        for (i = 0; i < 8; i++)
        {
            multivector_t mv = PGA3D_point(
                c.v[i].x,
                c.v[i].y,
                c.v[i].z
            );

            
            /*
            M = PGA3D_rotor(0.01, PGA3D_wedge(PGA3D_e2(), PGA3D_e3()));
            M = PGA3D_mul(M, PGA3D_rotor(0.01, PGA3D_wedge(PGA3D_e3(), PGA3D_e1())));
            */

            mv = PGA3D_mul(M, PGA3D_mul( mv, PGA3D_reverse(M)));

            c.v[i].x = mv.mvec[13];
            c.v[i].y = mv.mvec[12];
            c.v[i].z = mv.mvec[11];
        }
        
        for (i =0; i < 4; i++){
            SDL_RenderDrawLine(r,
                (c.v[i+4].x +1)* 400, (c.v[i+4].y +1)* 400,
                (c.v[(i+1)%4+4].x +1)* 400, (c.v[(i+1)%4+4].y +1)* 400
            );
        }       
        
        SDL_SetRenderDrawColor(r, 0xFF, 0xDD, 0x30, 0xff);
        for (i =0; i < 4; i++){
            SDL_RenderDrawLine(r,
                (c.v[i].x +1)* 400, (c.v[i].y +1)* 400,
                (c.v[i+4].x +1)* 400, (c.v[i+4].y +1)* 400
            );
        }       

        SDL_SetRenderDrawColor(r, 0xA0, 0x00, 0xFF, 0xff);

        SDL_RenderDrawLine(r,
            (c.v[0].x +1)* 400, (c.v[0].y +1)* 400,
            anchor_box.x + anchor_box.w/2,
            anchor_box.y + anchor_box.h/2
        );

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

multivector_t spring(multivector_t M, multivector_t B, multivector_t anchor_box, multivector_t corner, float k)
{
    multivector_t grav, hooke, damp;
    grav = PGA3D_dual(
        PGA3D_mul(
            PGA3D_mul(
                PGA3D_reverse(M),
                PGA3D_smul(
                    -0.981,
                    PGA3D_wedge(PGA3D_e0(), PGA3D_e2())
                )
            ),
            M
        )
    );

    hooke = PGA3D_smul(
        k,
        PGA3D_vee(
            PGA3D_mul(
                PGA3D_reverse(M),
                PGA3D_mul(
                    anchor_box,
                    M
                )
            ),
            corner
        )
    );

    damp = PGA3D_dual(PGA3D_smul(-0.25, B));

    return PGA3D_add(
        PGA3D_add(grav, hooke),
        damp
    );

}

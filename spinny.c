#include <stdio.h>
#include <string.h>
#include <math.h>
#include <SDL.h>
#include "pga3d.h"

/* Spring Force-Torque Function */
multivector_t spring_ftf(multivector_t M, multivector_t spring, multivector_t vertex, float k);
/* Gravity Force-Torque Function */
multivector_t gravity_ftf(multivector_t M, float grav_acc);
/* Damping Force-Torque Function */
multivector_t damping_ftf(multivector_t B, float damp_factor);

typedef struct vert{
    float x, y, z;
} vert;

typedef struct cube{
    vert v[8];
}cube;

int
render(
        SDL_Renderer *r, 
        cube * c,
        vert * a,
        vert lookat
);

vert 
viewport(
        SDL_Renderer *r, 
        vert lookat,
        float scale,
        vert wc
);

int main(void)
{
    SDL_Window *w = NULL;
    SDL_Renderer *r = NULL;
    SDL_Event  e;
    float dt, t, t_scale_factor;
    cube c;
    int i;
    vert a;
    multivector_t M, B, spring, dM, dB, net_ftf;

    c.v[0] = (vert) {-0.5, 0.5, 0.5};
    c.v[1] = (vert) {0.5, 0.5, 0.5};
    c.v[2] = (vert) {0.5, -0.5, 0.5};
    c.v[3] = (vert) {-0.5, -0.5, 0.5};
    
    c.v[4] = (vert) {-0.5, 0.5, -0.5};
    c.v[5] = (vert) {0.5, 0.5, -0.5};
    c.v[6] = (vert) {0.5, -0.5, -0.5};
    c.v[7] = (vert) {-0.5, -0.5, -0.5};
    
    a = (vert) {0, 0.5, 0};

    t_scale_factor = 1.0/100.0;

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
    memset(&B, 0, sizeof(B));
    spring = PGA3D_mul(
        //PGA3D_ssub(1, PGA3D_smul(0.5,PGA3D_wedge(PGA3D_e0(), PGA3D_e2()))),
        PGA3D_point(a.x, a.y, a.z),
        PGA3D_point(c.v[0].x,c.v[0].y,c.v[0].z)
    );
    

    do{
        SDL_SetRenderDrawColor(r, 0, 0, 0, 0xff);
        SDL_RenderClear(r);

        /* DEBUG: */
        /* 
        net_ftf = gravity_ftf(M, -9.81); 

        */
        net_ftf = PGA3D_add(
            gravity_ftf(M, -9.81),
            PGA3D_add(
                spring_ftf(M, spring, PGA3D_point(c.v[0].x,c.v[0].y,c.v[0].z), 12),
                damping_ftf(B, -0.25)
            )
        );


        dB = PGA3D_dual( 
            PGA3D_sub(
                net_ftf, 
                PGA3D_smul(0.5,
                    PGA3D_sub(
                        PGA3D_mul(PGA3D_dual(B), B),
                        PGA3D_mul(B, PGA3D_dual(B))
                    )
                )
            )
        );
        /*
        dB = PGA3D_smul(
            0.5,
            PGA3D_sub(
                PGA3D_mul(PGA3D_dual(B), B),
                PGA3D_mul(B, PGA3D_dual(B))
            )
        );
        */

        dM = PGA3D_mul(
            PGA3D_smul(-0.5, M),
            B
        );
        /* TODO: RK4 */
        /* Euler Integration */
        M = PGA3D_add(M, PGA3D_smul(t_scale_factor*dt/1000.0, dM));
        B = PGA3D_add(B, PGA3D_smul(t_scale_factor*dt/1000.0, dB));
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
        render(r, &c, &a, (vert){0,0,0}); 
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
int
render(
        SDL_Renderer *r, 
        cube * c,
        vert * a,
        vert lookat
)
{        
    int ww, wh, i;
    SDL_Rect anchor_box;
    vert pt0, pt1;
    SDL_GetRendererOutputSize(r, &ww, &wh);

    /* Set Anchor */
    anchor_box.x = ww/2*(a->x) - 4;
    anchor_box.y = wh/2*(a->y*-1) -4;
    anchor_box.w = 8;
    anchor_box.h = 8;

    /* Draw Anchor */
    SDL_SetRenderDrawColor(r, 0, 0xFF, 0x10, 0xff);
    SDL_RenderDrawRect(r, &anchor_box);

    /* Draw cube */
    SDL_SetRenderDrawColor(r, 0, 0xFF, 0x80, 0xff);

    for (i =0; i < 4; i++){
        pt0 = viewport(r, lookat, 1, c->v[i]);
        pt1 = viewport(r, lookat, 1, c->v[(i+4)%4]);
        SDL_RenderDrawLine(r,
            pt1.x,pt0.x,
            pt1.x,pt1.x
        );
    }

    for (i =0; i < 4; i++){
        pt0 = viewport(r, lookat, 1, c->v[i+4]);
        pt1 = viewport(r, lookat, 1, c->v[(i+4)%4+4]);
        SDL_RenderDrawLine(r,
            pt1.x,pt0.x,
            pt1.x,pt1.x
        );
    }       
    
    SDL_SetRenderDrawColor(r, 0xFF, 0xDD, 0x30, 0xff);
    for (i =0; i < 4; i++){
        pt0 = viewport(r, lookat, 1, c->v[i]);
        pt1 = viewport(r, lookat, 1, c->v[i+4]);
        SDL_RenderDrawLine(r,
            pt1.x,pt0.x,
            pt1.x,pt1.x
        );
    }       

    /* Draw spring */
    SDL_SetRenderDrawColor(r, 0xA0, 0x00, 0xFF, 0xff);

    pt0 = viewport(r, lookat, 1, c->v[0]);
    SDL_RenderDrawLine(r,
        pt0.x, pt1.y,
        anchor_box.x + anchor_box.w/2,
        anchor_box.y + anchor_box.h/2
    );

    SDL_RenderPresent(r);

    return 0;
}


multivector_t spring_ftf(multivector_t M, multivector_t spring, multivector_t vertex, float k)
{
    return PGA3D_smul(
        k,
        PGA3D_vee(
            PGA3D_mul(
                PGA3D_reverse(M),
                PGA3D_mul(
                    spring,
                    M
                )
            ),
            vertex
        )
    );
}

multivector_t gravity_ftf(multivector_t M, float grav_acc)
{
    return PGA3D_dual(
        PGA3D_mul(
            PGA3D_reverse(M),
            PGA3D_mul(
                PGA3D_smul(
                    grav_acc,
                    PGA3D_wedge(PGA3D_e0(), PGA3D_e2())
                ),
                M
            )
        )
    );
}

multivector_t damping_ftf(multivector_t B, float damp_factor)
{
    return PGA3D_dual(PGA3D_smul(damp_factor, B));
}



vert 
viewport(
        SDL_Renderer *r, 
        vert lookat,
        float scale,
        vert wc
){
    multivector_t M, pt, l;
    int ww, wh;
    SDL_GetRendererOutputSize(r, &ww, &wh);

    memset(&M, 0, sizeof(M));
    memset(&pt, 0, sizeof(pt));
    memset(&l, 0, sizeof(l));
    l = PGA3D_vee(
        PGA3D_point(ww/2, wh/2, 0),
        PGA3D_point(lookat.x, lookat.y, lookat.z)
    );
    M = PGA3D_mul(PGA3D_reverse(l),l);
    M = PGA3D_translator(sqrtf(M.mvec[0]), l);
    pt =  PGA3D_mul(
        M,
        PGA3D_mul(
            PGA3D_point(wc.x, wc.y, wc.z),
            PGA3D_reverse(M)
        )
    );
    return (vert){pt.mvec[13],pt.mvec[12],pt.mvec[11]};
}

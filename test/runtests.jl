using OrthoPoly
using Polynomials
using Test

function orthogonal_poly_test(F,Ps,num=[])
    if num == []
        num = ones(eltype(Ps[1]), length(Ps))
    end
    for n=0:length(Ps)-1
        f = F(n)
        p = Poly(Ps[n+1])/num[n+1]
        @test f == p
    end
end

println("Testing the first 11 Legendre polynomials")
orthogonal_poly_test(legendre, [[1],[0,1],[-1,0,3],[0,-3,0,5],[3,0,-30,0,35],[0,15,0,-70,0,63],[-5,0,105,0,-315,0,231],[0,-35,0,315,0,-693,0,429],[35,0,-1260,0,6930,0,-12012,0,6435],[0,315,0,-4620,0,18018,0,-25740,0,12155],[-63,0,3465,0,-30030,0,90090,0,-109395,0,46189]], [1 1 2 2 8 8 16 16 128 128 256])

println("Testing the first 16 Chebyshev polynomials of the first kind")
orthogonal_poly_test(chebyshevT, [[1],[0,1],[-1,0,2],[0,-3,0,4],[1,0,-8,0,8],[0,5,0,-20,0,16],[-1,0,18,0,-48,0,32],[0,-7,0,56,0,-112,0,64],[1,0,-32,0,160,0,-256,0,128],[0,9,0,-120,0,432,0,-576,0,256],[-1,0,50,0,-400,0,1120,0,-1280,0,512],[0,-11,0,220,0,-1232,0,2816,0,-2816,0,1024],[1,0,-72,0,840,0,-3584,0,6912,0,-6144,0,2048],[0,13,0,-364,0,2912,0,-9984,0,16640,0,-13312,0,4096],[-1,0,98,0,-1568,0,9408,0,-26880,0,39424,0,-28672,0,8192],[0,-15,0,560,0,-6048,0,28800,0,-70400,0,92160,0,-61440,0,16384]])

println("Testing the first 16 Chebyshev polynomials of the second kind")
orthogonal_poly_test(chebyshevU, [[1],[0,2],[-1,0,4],[0,-4,0,8],[1,0,-12,0,16],[0,6,0,-32,0,32],[-1,0,24,0,-80,0,64],[0,-8,0,80,0,-192,0,128],[1,0,-40,0,240,0,-448,0,256],[0,10,0,-160,0,672,0,-1024,0,512],[-1,0,60,0,-560,0,1792,0,-2304,0,1024],[0,-12,0,280,0,-1792,0,4608,0,-5120,0,2048],[1,0,-84,0,1120,0,-5376,0,11520,0,-11264,0,4096]])

println("Testing the first 11 Hermite polynomials")
orthogonal_poly_test(hermite_prob, [[1],[0,1],[-1,0,1],[0,-3,0,1],[3,0,-6,0,1],[0,15,0,-10,0,1],[-15,0,45,0,-15,0,1],[0,-105,0,105,0,-21,0,1],[105,0,-420,0,210,0,-28,0,1],[0,945,0,-1260,0,378,0,-36,0,1],[-945,0,4725,0,-3150,0,630,0,-45,0,1]])
orthogonal_poly_test(hermite, [[1],[0,2],[-2,0,4],[0,-12,0,8],[12,0,-48,0,16],[0,120,0,-160,0,32],[-120,0,720,0,-480,0,64],[0,-1680,0,3360,0,-1344,0,128],[1680,0,-13440,0,13440,0,-3584,0,256],[0,30240,0,-80640,0,48384,0,-9216,0,512],[-30240,0,302400,0,-403200,0,161280,0,-23040,0,1024]])

println("Testing the first 9 Laguerre polynomials")
orthogonal_poly_test(n -> laguerre(n, Poly([0,1//1])), [[1],[1,-1],[2,-4,1],[6,-18,9,-1],[24,-96,72,-16,1],[120,-600,600,-200,25,-1],[720,-4320,5400,-2400,450,-36,1],[5040,-35280,52920,-29400,7350,-882,49,-1],[40320,-322560,564480,-376320,117600,-18816,1568,-64,1]], [1 1 2 6 24 120 720 5040 40320]//1)

println("Testing the first 4 associated Laguerre polynomials")
for α=[1,1.0,2,2.0,2.1]
    orthogonal_poly_test(n -> laguerre_assoc(n,α,Poly([0,one(α)])),
                         [[1],[1+α,-1],[(α+2)*(α+1),-2(α+2),1],[(α+1)*(α+2)*(α+3),-3(α+2)*(α+3),3(α+3),-1]], [1.0 1 2 6])
end

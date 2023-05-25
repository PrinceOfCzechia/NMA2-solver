import numpy as np
# from matplotlib import pyplot as plt

# initial variables
N = 101 # 1000+ for a good fit
h = float( 1/N )
p = -1

X = np.linspace( 0, 1, N )

def q( x ):
    return 4 * np.pi**2 * x**2

Q = np.zeros( N )
for i in range( N ):
    Q[ i ] = q( X[i] )

def f( x ):
    return 2 * np.pi * np.cos( np.pi * x**2 )

F = np.zeros( N )
for i in range( N ):
    F[ i ] = f( X[i] )

# Thomas' algorithm to solve a tridiagonal linear system
def thomas( M: np.array, RHS: np.array ):
    forward_thomas( M, RHS )
    return backward_thomas( M1, RHS1 )

def forward_thomas( M: np.array, RHS: np.array ):
    global M1, RHS1
    M1 = M.copy()
    RHS1 = RHS.copy()
    for i in range( M1.shape[0] - 1 ):
        a = M1[ i ][ i ]
        b = M1[ i+1 ][ i ]
        M1[ i+1 ][ i ] = 0
        M1[ i+1 ][ i+1 ] -= b/a * M1[ i ][ i+1 ]
        RHS1[ i+1 ] -= b/a * RHS1[ i ]

def backward_thomas( M: np.array, RHS: np.array ):
    n = RHS.shape[0]
    u = np.zeros( n )
    u[ n-1 ] = RHS[ n-1 ] / M[ n-1 ][ n-1 ]
    for i in range( 1, n ):
        u[ n-i-1 ] = ( RHS[ n-i-1 ] - M[ n-i-1 ][ n-i ] * u[ n-i ] ) / M[ n-i-1][ n-i-1 ]
    return u

# setting up the linear system
def setup_matrix( N ):
    M = np.zeros( [ N, N ] )
    M[ 0 ][ 0 ] = 1
    for i in range( 1, N-1 ):
        M[ i ][ i-1 ] = -p / h**2
        M[ i ][ i ] = ( p+p ) / h**2 + Q[ i ]
        M[ i ][ i+1 ] = -p / h**2
    M[ N-1 ][ N-1 ] = 1
    return M

def setup_RHS( N, gamma_0, gamma_1 ):
    RHS = np.zeros( N )
    RHS[ 0 ] = gamma_0
    for i in range( 1, N-1 ):
        RHS[ i ] = F[ i ]
    RHS[ N-1 ] = gamma_1
    return RHS

# numerical solution
U = thomas( setup_matrix( N ), setup_RHS( N, 0, 0 ) )

# exact solution
PY = np.zeros( N )
for i in range( N ):
    PY[ i ] = np.sin( np.pi * X[ i ]**2)

# comparison
def l2_norm( X: np.array):
    sum = 0.0
    for i in range( U.shape[ 0 ] ):
        sum += X[ i ]**2
    return np.sqrt( sum )

# error = l2_norm(PY-U)
# avg_error = error/N
# print( 'Error =', error, 'Average error =', avg_error )

# print gnuplot-compatible result
f = open( 'solution.txt', 'a' )
for i in range( N ):
    # print( X[ i ], ' ', U[ i ] )
    f.write( str(X[ i ]) +  '  ' + str(U[ i ]) + '\n' )
    
g = open( 'exact.txt', 'a' )
for i in range( N ):
    g.write( str(X[ i ]) +  '  ' + str(PY[ i ]) + '\n' )


# uncomment for python plots
# plt.plot( X, PY, 'b-', label = 'exact solution' )
# plt.plot( X, U, '--', color = 'orange', label = 'numerical solution' )
# plt.show()
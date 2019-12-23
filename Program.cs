// Based on http://cliffle.com/p/dangerust/1/
// Which is a direct port of https://benchmarksgame-team.pages.debian.net/benchmarksgame/program/nbody-gcc-8.html
// Takes about 6.3 seconds on my computer (and the Rust solution takes 3.3)

using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

[StructLayout(LayoutKind.Sequential)]
struct Vector3
{
    public double X;
    public double Y;
    public double Z;

    public Vector3(double x, double y, double z)
    {
        X = x;
        Y = y;
        Z = z;
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3 operator +(in Vector3 a, in Vector3 b) {
        return new Vector3(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3 operator -(in Vector3 a, in Vector3 b) {
        return new Vector3(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3 operator *(in Vector3 a, in Vector3 b) {
        return new Vector3(a.X * b.X, a.Y * b.Y, a.Z * b.Z);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3 operator *(in Vector3 a, double b) {
        return new Vector3(a.X * b, a.Y * b, a.Z * b);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3 operator *(double a, in Vector3 b) {
        return new Vector3(a * b.X, a * b.Y, a * b.Z);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double SumOfSquares() {
        return X * X + Y * Y + Z * Z;
    }

    public override string ToString()
    {
        return $"Vector3({X}, {Y}, {Z})";
    }
}

struct Body
{
    public Vector3 Position;
    public Vector3 Velocity;
    public double Mass;
}

public static class NBody
{
    const double SOLAR_MASS = 4.0 * Math.PI * Math.PI;
    const double DAYS_PER_YEAR = 365.24;
    const int BODIES_COUNT = 5;
    const int INTERACTIONS_COUNT = BODIES_COUNT * (BODIES_COUNT - 1) / 2;
    const int ROUNDED_INTERACTIONS_COUNT = INTERACTIONS_COUNT + INTERACTIONS_COUNT % 2;

    public static void Main(String[] args)
    {
        int n = args.Length > 0 ? Int32.Parse(args[0]) : 5000;
        var solarBodies = new Body[] {
            // Sun,
            new Body { Mass = SOLAR_MASS },
            // Jupiter
            new Body {
                Position = new Vector3(
                    4.84143144246472090e+00,
                    -1.16032004402742839e+00,
                    -1.03622044471123109e-01
                ),
                Velocity = new Vector3(
                    1.66007664274403694e-03 * DAYS_PER_YEAR,
                    7.69901118419740425e-03 * DAYS_PER_YEAR,
                    -6.90460016972063023e-05 * DAYS_PER_YEAR
                ),
                Mass = 9.54791938424326609e-04 * SOLAR_MASS
            },
            // Saturn
            new Body {
                Position = new Vector3(
                    8.34336671824457987e+00,
                    4.12479856412430479e+00,
                    -4.03523417114321381e-01
                ),
                Velocity = new Vector3(
                    -2.76742510726862411e-03 * DAYS_PER_YEAR,
                    4.99852801234917238e-03 * DAYS_PER_YEAR,
                    2.30417297573763929e-05 * DAYS_PER_YEAR
                ),
                Mass = 2.85885980666130812e-04 * SOLAR_MASS
            },
            // Uranus
            new Body {
                Position = new Vector3(
                    1.28943695621391310e+01,
                    -1.51111514016986312e+01,
                    -2.23307578892655734e-01
                ),
                Velocity = new Vector3(
                    2.96460137564761618e-03 * DAYS_PER_YEAR,
                    2.37847173959480950e-03 * DAYS_PER_YEAR,
                    -2.96589568540237556e-05 * DAYS_PER_YEAR
                ),
                Mass = 4.36624404335156298e-05 * SOLAR_MASS
            },
            // Neptune
            new Body {
                Position = new Vector3(
                    1.53796971148509165e+01,
                    -2.59193146099879641e+01,
                    1.79258772950371181e-01
                ),
                Velocity = new Vector3(
                    2.68067772490389322e-03 * DAYS_PER_YEAR,
                    1.62824170038242295e-03 * DAYS_PER_YEAR,
                    -9.51592254519715870e-05 * DAYS_PER_YEAR
                ),
                Mass = 5.15138902046611451e-05 * SOLAR_MASS
            },
        };

        OffsetMomentum(solarBodies);
        Console.WriteLine("{0:F9}", CalculateEnergy(solarBodies));

        for (var i = 0; i < n; i++) {
            Advance(solarBodies);
        }

        Console.WriteLine("{0:F9}", CalculateEnergy(solarBodies));
    }

    private static void OffsetMomentum(Body[] bodies)
    {
        ref var sun = ref bodies[0];

        for (var i = 0; i < bodies.Length; i++) {
            ref var body = ref bodies[i];
            var relativeMass = body.Mass / SOLAR_MASS;
            sun.Velocity -= body.Velocity * relativeMass;
        }
    }

    private static double CalculateEnergy(Body[] bodies) {
        var energy = 0.0;

        for (var i = 0; i < bodies.Length; i++) {
            ref var body = ref bodies[i];
            energy += 0.5 * body.Mass * body.Velocity.SumOfSquares();

            for (var j = i + 1; j < bodies.Length; j++) {
                ref var otherBody = ref bodies[j];
                var positionDelta = body.Position - otherBody.Position;
                energy -= (body.Mass * otherBody.Mass) / Math.Sqrt(positionDelta.SumOfSquares());
            }
        }

        return energy;
    }

    [StructLayout(LayoutKind.Sequential, Pack = 16)]
    unsafe struct Align16 {
        public fixed double Data[ROUNDED_INTERACTIONS_COUNT]; 
    }

    private static unsafe void Advance(Body[] bodies) {
        var positionDeltas = stackalloc Align16[3];
        Align16 magnitudes;

        var k = 0;

        for (var i = 0; i < bodies.Length - 1; i++) {
            ref var body = ref bodies[i];

            for (var j = i + 1; j < bodies.Length; j++) {
                ref var otherBody = ref bodies[j];
                positionDeltas[0].Data[k] = body.Position.X - otherBody.Position.X;
                positionDeltas[1].Data[k] = body.Position.Y - otherBody.Position.Y;
                positionDeltas[2].Data[k] = body.Position.Z - otherBody.Position.Z;
                k++;
            }
        }

        for (var i = 0; i < ROUNDED_INTERACTIONS_COUNT / 2; i++) {
            var positionDelta = stackalloc Vector128<double>[3];
            positionDelta[0] = ((Vector128<double>*)(positionDeltas + 0))[i];
            positionDelta[1] = ((Vector128<double>*)(positionDeltas + 1))[i];
            positionDelta[2] = ((Vector128<double>*)(positionDeltas + 2))[i];

            var distanceSquared = Sse2.Add(
                Sse2.Add(
                    Sse2.Multiply(positionDelta[0], positionDelta[0]),
                    Sse2.Multiply(positionDelta[1], positionDelta[1])
                ),
                Sse2.Multiply(positionDelta[2], positionDelta[2])
            );

            var distanceReciprocal = Sse2.ConvertToVector128Double(
                Sse.ReciprocalSqrt(Sse2.ConvertToVector128Single(distanceSquared))
            );

            for (var x = 0; x < 2; x++) {
                distanceReciprocal = Sse2.Subtract(
                    Sse2.Multiply(distanceReciprocal, Vector128.Create(1.5, 1.5)),
                    Sse2.Multiply(
                        Sse2.Multiply(
                            Sse2.Multiply(
                                Vector128.Create(0.5, 0.5),
                                distanceSquared
                            ),
                            distanceReciprocal
                        ),
                        Sse2.Multiply(
                            distanceReciprocal,
                            distanceReciprocal
                        )
                    )
                );
            }

            ((Vector128<double>*)&magnitudes)[i] = Sse2.Multiply(
                Sse2.Divide(Vector128.Create(0.01), distanceSquared),
                distanceReciprocal
            );
        }

        k = 0;
        for (var i = 0; i < bodies.Length - 1; i++) {
            ref var body = ref bodies[i];
            for (var j = i + 1; j < bodies.Length; j++) {
                ref var otherBody = ref bodies[j];

                var iMassMagnitude = body.Mass * magnitudes.Data[k];
                var jMassMagnitude = otherBody.Mass * magnitudes.Data[k];

                var a = positionDeltas[0].Data[k];
                var b = positionDeltas[1].Data[k];
                var c = positionDeltas[2].Data[k];

                body.Velocity.X -= a * jMassMagnitude;
                otherBody.Velocity.X += a * iMassMagnitude;

                body.Velocity.Y -= b * jMassMagnitude;
                otherBody.Velocity.Y += b * iMassMagnitude;

                body.Velocity.Z -= c * jMassMagnitude;
                otherBody.Velocity.Z += c * iMassMagnitude;
                k++;
            }
        }

        for (var i = 0; i < bodies.Length; i++) {
            bodies[i].Position += 0.01 * bodies[i].Velocity;
        }
    }
}

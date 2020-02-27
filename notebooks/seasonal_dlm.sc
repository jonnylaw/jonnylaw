import $ivy.`org.scalanlp::breeze:0.13.2`
import breeze.linalg._
import breeze.stats.distributions._
import breeze.numerics._
import java.io.{File, PrintWriter}

case class Data(time: Double, observation: Double, state: DenseVector[Double]) {
  override def toString = s"$time, $observation, ${state.data.mkString(", ")}"
}
case class Dlm(f: DenseMatrix[Double], g: DenseMatrix[Double])

/**
* Simulate a single step from a dlm
*/
def simStep(p: Parameters, model: Dlm): Data => Rand[Data] = d => {
  for {
    w <- MultivariateGaussian(DenseVector.zeros(p.w.cols), p.w)
    x1 = model.g * d.state + w
    v <- Gaussian(0, p.v)
    y = model.f.toDenseVector dot x1 + v
  } yield Data(d.time + 1.0, y, x1)
}

/**
  * Simulate from the DLM
  * @param p the parameters of the DLM
  * @param init_state the initial state of the DLM
  * @return a Process object representing the DLM as a MarkovChain
  */
def simMarkov(p: Parameters, model: Dlm): Rand[Process[Data]] = {
  for {
    x0 <- MultivariateGaussian(p.m0, p.c0)
    y0 <- Gaussian(model.f.toDenseVector dot x0, p.v)
    init = Data(0.0, y0, x0)
  } yield MarkovChain(init)(simStep(p, model))
}

/**
  * The free parameters of a DLM, with constant observation and system variance
*/
case class Parameters(
  v: Double,
  w: DenseMatrix[Double],
  m0: DenseVector[Double],
  c0: DenseMatrix[Double])

/**
  * Block concatenate two matrices into a larger matrix
  */
def blockConcat(x: DenseMatrix[Double], y: DenseMatrix[Double]): DenseMatrix[Double] = {
  val res = DenseMatrix.zeros[Double](x.rows + y.rows, x.cols + y.cols)
  res(0 to (x.rows - 1), 0 to (x.cols - 1)) := x
  res(x.rows to (x.rows + y.rows - 1), x.cols to (x.cols + y.cols - 1)) := y
  res
}

/**
  * Create a 2x2 rotation matrix with frequency = 2 * pi / period
  * @param period, the period of the rotation
  */
def rotationMatrix(period: Int, harmonic: Int): DenseMatrix[Double] = {
  val freq = 2 * math.Pi * harmonic / period
  DenseMatrix((cos(freq), -sin(freq)), (sin(freq), cos(freq)))
}

/**
  * A fourier-form seasonal DLM
  */
def SeasonalDlm(period: Int, harmonics: Int): Dlm = {
  def f = DenseMatrix.tabulate(1, harmonics * 2)((i, j) => if (j % 2 == 0) 1.0 else 0.0)

  def g = (1 to harmonics).
  map (h => rotationMatrix(period, h)).
  reduce((a, b) => blockConcat(a, b))

  Dlm(f, g)
}

val p = Parameters(
  3.0,
  DenseMatrix.eye[Double](6),
  DenseVector.zeros[Double](6),
  DenseMatrix.eye[Double](6)*10.0)

val model = SeasonalDlm(24, 3)

val sims = simMarkov(p, model).
  draw.
  steps.
  take(24*7)

val pw = new PrintWriter("notebooks/data/SeasonalDlm.csv")
pw.write(sims.mkString("\n"))
pw.close()

import $ivy.`org.scalanlp::breeze:0.13.2`
import breeze.stats.distributions._
import breeze.linalg._
import breeze.numerics.exp
import java.io.{PrintWriter, File}

case class Parameters(mu: DenseVector[Double], sigma: Double)

def model(params: Parameters) =
    MultivariateGaussian(params.mu, diag(DenseVector.fill(2)(params.sigma)))

// Simulate values from the model.

val p = Parameters(DenseVector(2.0, 3.0), 0.5)
val xy = model(p).sample(100)

val pw1 = new PrintWriter(new File("notebooks/data/BivariateSimulated.csv"))
    xy.
      foreach( d =>
        pw1.write(s"${d(0)}, ${d(1)}\n")
      )
    pw1.close()

// Specify the likelihood and prior distributions and define the log measure.

def likelihood(points: Seq[DenseVector[Double]], p: Parameters) =
    points.map { point =>
      MultivariateGaussian(p.mu, diag(DenseVector.fill(2)(p.sigma))).logPdf(point)
    }.reduce((x, y) => x + y)

def prior(p: Parameters) = {
    MultivariateGaussian(DenseVector(2.0, 3.0), diag(DenseVector.fill(2)(3.0))).logPdf(p.mu) +
      Gamma(shape = 0.5, scale = 2.0).logPdf(1/p.sigma)
}

val logMeasure = (p: Parameters) => likelihood(xy, p) + prior(p)

def simPrior = for {
    mu <- MultivariateGaussian(DenseVector(2.0, 3.0), diag(DenseVector.fill(2)(3.0)))
    sigma <- Gamma(2.0, 3.0)
} yield Parameters(mu, sigma)

def propose(scale: Double)(p: Parameters) =
    for {
      innov <- MultivariateGaussian(DenseVector.fill(3)(0.0), diag(DenseVector.fill(3)(scale)))
      mu = p.mu + innov(0 to 1)
      sigma = p.sigma * exp(innov(2))
    } yield Parameters(mu, sigma)

val params = MarkovChain.metropolis(simPrior.draw, propose(0.025))(logMeasure).
    steps.
    take(10000).
    toSeq.
    map(p => s"${p.mu(0)}, ${p.mu(1)}, ${p.sigma}")

val pw = new PrintWriter(new File("notebooks/data/Parameters.csv"))
    params.
      take(10000).
      foreach { p =>
        pw.write(p + "\n")
      }
    pw.close()

import $ivy.`com.stripe::rainier-core:0.2.2`
import com.stripe.rainier.core._
import com.stripe.rainier.sampler._
import java.io.{File, PrintWriter}

val (alpha, beta, sigma) =  (-1.5, 2.0, 0.5)

// using sampling monad to create synthetic data
val lm = for {
    x <- Normal(0, 1).param
    y <- Normal(alpha + beta * x, sigma).param
} yield (x, y)

implicit val s = RNG.default
val sims = lm.sample(100)

val pw = new PrintWriter(new File("notebooks/data/lm_sims.csv"))
pw.write("x, y\n")
pw.write(sims.map { case (x, y) => s"$x, $y" }.mkString("\n"))
pw.close()

import com.stripe.rainier.compute._
// import com.stripe.rainier.core._
// import com.stripe.rainier.sampler._

def linearModel(data: Seq[(Double, Double)]): RandomVariable[Map[String, Real]] = for {
    alpha <- Normal(0, 5).param
    beta <- Normal(0, 5).param
    sigma <- LogNormal(2, 2).param
    _ <- Predictor[Double].from { x =>
      Normal(alpha + beta * x, sigma)
    }
    .fit(data)
  } yield Map("alpha" -> alpha, "beta" -> beta, "sigma" -> sigma)

val iters = linearModel(sims).sample(HMC(5), 5000, 100000, 100)


val pw1 = new PrintWriter(new File("notebooks/data/lm_params.csv"))
pw1.write(s"${iters.head.keys.mkString(", ")}\n")
pw1.write(iters.
  map(ps => ps.values.mkString(", ")).mkString("\n"))
pw1.close()

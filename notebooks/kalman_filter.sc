import $ivy.`org.scalanlp::breeze:0.13.2`
import breeze.stats.distributions._
import java.io.{PrintWriter, File}

case class Data(time: Double, observation: Double, state: Option[Double]) {
  override def toString = s"${time}, ${observation}, ${state.get}"
}
case class Parameters(v: Double, w: Double, m0: Double, c0: Double)

def stepRw(p: Parameters): Double => Rand[Double] =
  x => Gaussian(x, p.w)

def simulate(p: Parameters, n: Int): Seq[Data] = {
  val x0 = Gaussian(p.m0, math.sqrt(p.c0)).draw
  val stateSpace = Seq.iterate(x0, n)(x => stepRw(p)(x).draw)
  stateSpace.zipWithIndex map { case (x, t) =>
    Data(t, x + Gaussian(0, math.sqrt(p.v)).draw, Some(x))
  }
}

val p = Parameters(3.0, 0.5, 0.0, 10.0)

// simulate 16 different realisations of 100 observations, representing 16 stations
val data: Map[Int, Seq[Data]] = (1 to 16).
  map((id: Int) => (id, simulate(p, 100))).
  toMap

val pw = new PrintWriter(new File("notebooks/data/firstOrderDlm.csv"))
pw.write(
  data.
    flatMap{ case (id, data) =>
      data map (x => id + ", " + x.toString)}.
    mkString("\n"))
pw.close()

// state of the filter containing the data, parameters and log-likelihood
case class FilterState(data: Data, p: Parameters, ll: Double) {
  override def toString = s"${data.toString}, ${p.m0}, ${p.c0}"
}

def filter(p: Parameters): (FilterState, Data) => FilterState = (s, d) => {
  val r = s.p.c0 + p.w
  val q = r + p.v
  val e = d.observation - s.p.m0

  // kalman gain
  val k = r / q
  val c1 = k * p.v
  val m1 = s.p.m0 + k * e
  val ll = Gaussian(m1, math.sqrt(c1)).logPdf(d.observation)

  // return the data with the expectation of the hidden state and the updated Parameters
  FilterState(Data(d.time, d.observation, Some(m1)), Parameters(p.v, p.w, m1, c1), s.ll + ll)
}

def filterSeries(data: Seq[Data], initParams: Parameters): Seq[FilterState] = {
  val initFilter = FilterState(data.head, initParams, 0.0) // initialise the filter

  data.
    scanLeft(initFilter)(filter(initParams)).
    drop(1)
}

// Now, we can apply the filter to all the stations simultaneously.

val filtered = data.
  mapValues(filterSeries(_, p)) // apply the filter to the sorted data

// write the filter for all stations to a file
val pw1 = new PrintWriter(new File("notebooks/data/filteredDlm.csv"))
pw1.write(filtered.
           flatMap{ case (id, data) =>
             data map (x => id + ", " + x.toString)}.
           mkString("\n"))
pw1.close()

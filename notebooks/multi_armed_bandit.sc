import $ivy.`org.scalanlp::breeze:0.13.2`

import breeze.stats.distributions._
import breeze.linalg._
import java.io.{File, PrintWriter}

case class BanditState(
  reward: Array[Double],
  longTermReward: List[Double],
  actions: Map[Int, Int]
)

def banditStep(
  epsilon: Double,
  reward: Int => Rand[Double])(s: BanditState): Rand[BanditState] = {
  for {
    u <- Uniform(0, 1)
    nextAction <- if (u < epsilon) {
      Multinomial(DenseVector.ones[Double](s.actions.size))
    } else {
      Rand.always(s.longTermReward.zipWithIndex.maxBy(_._1)._2)
    }
    newReward <- reward(nextAction)
    prevCount = s.actions.get(nextAction).get
    nextCount  = prevCount + 1
    newLongTermReward = s.longTermReward(nextAction) + (newReward - s.longTermReward(nextAction)) / nextCount
  } yield BanditState(s.reward :+ newReward,
                      s.longTermReward.updated(nextAction, newLongTermReward),
                      s.actions.updated(nextAction, nextCount))
}

def buildActions(actions: Int): Map[Int, Int] = {
  (0 until actions).map(a => a -> 0).toMap
}

def bandit(
  epsilon: Double,
  actions: Int,
  reward: Int => Rand[Double],
  n: Int): BanditState = {

  val initState = BanditState(Array(0.0), List.fill(10)(0.0), buildActions(actions))
  MarkovChain(initState)(banditStep(epsilon, reward)).steps.drop(n-1).next
}

val qs = Gaussian(0, 1).sample(10)

// The reward is selected from a N(q(A_t), 1)
def r(qa: Seq[Double])(action: Int): Rand[Double] =
  Gaussian(qa(action), 1)

// plot distribution of rewards
val oneBandit = bandit(0.5, 10, r(qs), 1000)
val pw1 = new PrintWriter(new File("notebooks/data/action_distribution.csv"))
pw1.write(oneBandit.actions.
            map { case (action, count) => s"$action, $count" }.
            mkString("\n"))
pw1.close()



// Calculate the average reward for n 10-arm bandit models with steps and an epsilon-greedy method
def averageReward(n: Int, steps: Int, epsilon: Double) = {
  Vector.fill(n)(DenseVector(bandit(epsilon, 10, r(qs), steps).reward)).
    reduce(_ + _).map(_ / n)
}

// average reward for eps = 0.1
val aveReward = averageReward(2000, 1000, 0.1)

// calculate the average reward for exp = 0.0, 0.1, 0.5
val eps = List(0.0, 0.1, 0.5)
val data = eps.map( eps => averageReward(2000, 1000, eps))

val pw = new PrintWriter(new File("notebooks/data/average_reward.csv"))
pw.write("epsilon, step, reward\n")
pw.write(data.
           zip(eps).
           flatMap { case (rs, e) => (List.fill(2000)(e), List.range(0, 2000), rs.data).
                        zipped.
                        map { case (es, t, r) => s"$es, $t, $r"} }.
           mkString("\n"))
pw.close()

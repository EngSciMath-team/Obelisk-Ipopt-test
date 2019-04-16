import org.coinor.Ipopt
import math._

case class Resource(id: Int, name: String, measurementUnit: String, naturalProduction: Double)

case class ResourceProduction(resource: Resource, production: Double)

case class Recipe(id: Int, name: String, production: List[ResourceProduction], utility: Double)

case class RecipeResourceProduction(resourceProduction: ResourceProduction, recipeId: Int)

case class RecipeUtility(recipe: Recipe, utility: Double)

case class RecipeSolution(recipeName: String, solution: Double)

case class SolverResult(objectiveValue: Double, recipeSolutions: Seq[RecipeSolution])

/* Extend the "Ipopt" class given by the Java interface. 
   This includes defining functions declared in Ipopt which modify certian input variables, 
   as opposed to returning the result.
   The error I get on trying to run this is that the function inputs aren't declared. */

class PlanIpopt(rec: Seq[Recipe], res: Array[Resource]) extends Ipopt{

  val recipes = rec
  val resources = res

  val index_style : Int = Ipopt.C_STYLE
  val numRecipes = recipes.length
  val numResources = resources.length
  this.create(numRecipes, numResources, numRecipes, numRecipes, index_style)

  // Set bounds on intensities (x), and the total surplus amount of a resource being produced (g).

  def get_bounds_info(n: Int, x_L: Array[Double], x_U: Array[Double],
                      m: Int, g_L: Array[Double], g_U: Array[Double]): Boolean {
    x_L = Array.fill(n)(0)
    g_L = Array.fill(m)(0)
    x_U = Array.fill(n)(2e19) 
    g_U = this.resources.map(resource => resource.naturalProduction)
    true
  }

  /* Give the starting intensities as (1,1,1,1,1...) for now.
     This will change once old intensities are passed to solver. */

  def get_starting_point(n: Int, init_x: Boolean, x: Array[Double],
                         init_z: Boolean, z_L: Array[Double], z_U: Array[Double],
                         m: Int, init_lambda: Boolean, lambda: Array[Double]): Boolean {
    require(init_z == false)
    require(init_lambda == false)
    x = Array.fill(n)(1)
    true
  }

  /* Define the function which gives the total utility of a set of intensities. */

  def eval_f(n: Int, x: Array[Double], new_x: Boolean, obj_value: Array[Double]) 
            : Boolean {
    obj_value = x.zip(this.recipes.map(_.utility).toArray).map(thisX => if (thisX._2 == 0) {
      0
    } else {
      thisX._2*log(thisX._1)
    }).sum
    true
  }

  /* Define the function which gives the *marginal* utilities of a set of intensities. 
     That is, the gradiant of the total utility */

  def eval_grad_f(n: Int, x: Array[Double], new_x: Boolean, grad_f: Array[Double]) 
                 : Boolean {
    grad_f = x.zip(this.recipes.map(_.utility).toArray).map(thisX => if (thisX._2 == 0) {
      0
    } else {
      thisX._2/thisX._1
    })
    true
  }

  /* Define the function which gives the surplus amount of each resource produced. */

  def eval_g(n: Int, x: Array[Double], new_x: Boolean, m: Int, g: Array[Double])
            : Boolean {
    this.recipes.foreach(recipe =>
      recipe.production.foreach(resourceProduction =>
        resId = resourceProduction.resource.id-1
        g(resId) += resourceProduction.production*x(resId)
    )
    true
  }

  /* Define the Jacobian matrix of the constraints.
     This corresponds in our case to the constraint matrix.
     iRow(n) and jCol(n) correspond to the coordinates of the nth entry */

  def eval_jac_g(n: Int, x: Array[Double], new_x: Boolean, m: Int, nele_jac: Int,
                             iRow: Array[Int], jCol: Array[Int], values Array[Double])
            : Boolean {

    var elementCounter : Int = 0

    this.recipes.foreach(recipe =>
      iRow(elementCounter) = recipe.id
      recipe.production.foreach(resourceProduction =>
        resId = resourceProduction.resource.id
        jCol(elementCounter) = resId
        values(resId) += resourceProduction.production
        elementCounter += 1
    )
    true
  }

  /* Define the Hessian matrix of the Lagrangian form of the problem.
     Not sure how to compactly explain this. */

  def eval_h(n: Int, x: Array[Double], new_x: Boolean, obj_factor: Double,
                         m: Int, lambda: Array[Double], new_lambda: Boolean,
                         nele_hess: Int, iRow: Array[Int], jCol: Array[Int], values: Array[Double])
            : Boolean {
    iRow = (0 to (m-1)).toArray
    jCol = iRow

    this.recipes.foreach(recipe =>
      recipe.production.foreach(resourceProduction =>
        resId = resourceProduction.resource.id
        values(resId) += -resourceProduction.production/pow(x(resId), 2)
    )
    true
  }
}

class Solver {

  def solve(recipes: Seq[Recipe]): SolverResult = {

    // Figure out how many resources there are in the data. If the number is readily available elsewhere, make it an input of solve and remove this.

    var nResources = 0

    recipes.foreach(recipe => 
      recipe.production.foreach(resourceProduction =>
        if (resourceProduction.resource.id > nResources) (
          nResources = resourceProduction.resource.id
        )
      )
    )

    // Process the incoming recipe case classes into an array of resource case classes.

    var resources : Array[Resource] = Array.ofDim(nResources)

    recipes.foreach(recipe =>
      recipe.production.foreach(resourceProduction =>
        resources(resourceProduction.resource.id-1) = resourceProduction.resource
      )
    )

    // Use the PlanIpopt class.

    PlanIpopt plan = new PlanIpopt(recipes, resources)
    x = Array.fill(recipes.length)(1)
    println("Status" ++ plan.OptimizeNLP(x))

    SolverResult(
      objectiveValue = plan.getObjectiveValue,
      recipeSolutions = (plan.getVariableValues.toList, recipes.map(_.name)).zipped.map(
        p => RecipeSolution(recipeName = p._2, solution = p._1)
      )
    )
  }
}


  def testDataOne : List[Recipe] = {
    // Resources
    val waterResource = Resource(id = 1, name = "Water", measurementUnit = "Cup", naturalProduction = 1)
    val iceResource = Resource(id = 2, name = "Ice", measurementUnit = "Cube", naturalProduction = 0)
    val potTimeResource = Resource(id = 3, name = "Pot Time", measurementUnit = "Pot Month", naturalProduction = 1)
    val flowerResource = Resource(id = 4, name = "Flower", measurementUnit = "Item", naturalProduction = 0)

    // Recipes
    val freezingRecipe = Recipe(
      id = 1,
      name = "Freezing",
      production = List(
        ResourceProduction(
          resource = waterResource,
          production = -2,
        ),
        ResourceProduction(
          resource = iceResource,
          production = 3,
        )
      ),
      utility = 0
    )

    val iceConsumptionRecipe = Recipe(
      id = 2,
      name = "Ice Consumption",
      production = List(
        ResourceProduction(
          resource = iceResource,
          production = -1,
        )
      ),
      utility = 1
    )

    val flowerGrowingRecipe = Recipe(
      id = 3,
      name = "Flower Growing",
      production = List(
        ResourceProduction(
          resource = waterResource,
          production = -1,
        ),
        ResourceProduction(
          resource = potTimeResource,
          production = -3,
        ),
        ResourceProduction(
          resource = flowerResource,
          production = 1,
        ),
      ),
      utility = 0
    )

    val flowerConsumptionRecipe = Recipe(
      id = 4,
      name = "Flower Consumption",
      production = List(
        ResourceProduction(
          resource = flowerResource,
          production = -1,
        )
      ),
      utility = 2
    )

  List(freezingRecipe,iceConsumptionRecipe, flowerGrowingRecipe, flowerConsumptionRecipe)

  }

  def testDataTwo : List[Recipe] = {
    // Resources
    val waterResource = Resource(id = 1, name = "Water", measurementUnit = "Cup", naturalProduction = 1)
    val iceResource = Resource(id = 2, name = "Ice", measurementUnit = "Cube", naturalProduction = 0)
    val potTimeResource = Resource(id = 3, name = "Pot Time", measurementUnit = "Pot Month", naturalProduction = 1)
    val flowerResource = Resource(id = 4, name = "Flower", measurementUnit = "Item", naturalProduction = 0)
    val carbonDioxideResource = Resource(id = 5, name = "Carbon Dioxide", measurementUnit="Gram", naturalProduction = 300)

    // Recipes
    val freezingRecipe = Recipe(
      id = 1,
      name = "Freezing",
      production = List(
        ResourceProduction(
          resource = waterResource,
          production = -2,
        ),
        ResourceProduction(
          resource = iceResource,
          production = 3,
        )
      ),
      utility = 0
    )

    val iceConsumptionRecipe = Recipe(
      id = 2,
      name = "Ice Consumption",
      production = List(
        ResourceProduction(
          resource = iceResource,
          production = -1,
        )
      ),
      utility = 1
    )

    val flowerGrowingRecipe = Recipe(
      id = 3,
      name = "Flower Growing",
      production = List(
        ResourceProduction(
          resource = waterResource,
          production = -1,
        ),
        ResourceProduction(
          resource = potTimeResource,
          production = -3,
        ),
        ResourceProduction(
          resource = flowerResource,
          production = 1,
        ),
        ResourceProduction(
          resource = carbonDioxideResource,
          production = -100
        )
      ),
      utility = 0
    )

    val flowerConsumptionRecipe = Recipe(
      id = 4,
      name = "Flower Consumption",
      production = List(
        ResourceProduction(
          resource = flowerResource,
          production = -1,
        )
      ),
      utility = 2
    )

  List(freezingRecipe,iceConsumptionRecipe, flowerGrowingRecipe, flowerConsumptionRecipe)

  }
}

object Main extends App {

  val solver = new Solver()
  val result = solver.solve(solver.testDataOne)

  println(result.toString)
}

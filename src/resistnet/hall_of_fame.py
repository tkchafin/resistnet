import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings

warnings.simplefilter("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
pd.options.mode.chained_assignment = None


class HallOfFame:
    """
    Class to represent a 'Hall of Fame' for storing and managing top-performing
    models.

    This class is designed to keep track of the best models based on various
    criteria such as fitness, log-likelihood, R-squared, AIC, etc. It allows
    for the addition and comparison of new models against existing ones.

    Attributes:
        data (pd.DataFrame): A DataFrame to store model attributes.
        variables (list): List of variable names used in the models.
        max_size (int): Maximum number of models to be stored in the hall of
                        fame.
        min_fitness (float): The minimum fitness value of the stored models.
        rvi (type): Description of rvi attribute.
        maw (type): Description of maw attribute.
        best (type): Description of best attribute.
        zero_threshold (float): Threshold value to consider as zero.
    """

    def __init__(self, variables, max_size, init_pop=None):
        """
        Initialize the Hall of Fame class with variables, maximum size, and an
        initial population.

        Args:
            variables (list): The list of variables to include in the model.
            max_size (int): The maximum size of the hall of fame.
            init_pop (list, optional): Initial population of models to consider
                                       Defaults to None.
        """
        cols = ["fitness"]
        for v in variables:
            cols.extend([str(v), f"{v}_weight", f"{v}_trans", f"{v}_shape", f"{v}_asym"])
        cols.extend(["loglik", "r2m", "aic", "delta_aic_null"])
        self.data = pd.DataFrame(columns=cols)
        self.variables = variables
        self.max_size = int(max_size)
        self.min_fitness = float("-inf")
        self.rvi = None
        self.maw = None
        self.best = None
        self.zero_threshold = 1e-17

        if init_pop is not None:
            self.check_population(init_pop)

    def check_population(self, pop):
        """
        Check and update the hall of fame with a new population.

        This method updates the hall of fame with models from the new
        population if they have better fitness scores than the existing ones.
        It also ensures that the size of the hall of fame does not exceed the
        maximum limit.

        Args:
            pop (list): A list of models to be considered for inclusion in the
                        hall of fame.
        """
        # only consider models w/ non neg-inf fitnesses
        popDF = pd.DataFrame(pop, columns=self.data.columns)
        popDF = popDF[popDF.fitness > float("-inf")]
        if popDF.shape[0] < 1:
            return

        popDF = popDF.sort_values("fitness", ascending=False)
        popDF = popDF.drop_duplicates(keep="first", ignore_index=True)
        popDF = popDF.reset_index(drop=True)
        space = self.max_size - self.data.shape[0]

        # Drop columns with all NA values before concatenation
        popDF = popDF.dropna(axis=1, how='all')

        if space > 0:
            select_size = min(space, popDF.shape[0])
            self.data = pd.concat(
                [self.data, popDF[:select_size]], ignore_index=True
            )
            self._update_data_frame()
        else:
            if popDF["fitness"].max() > self.min_fitness:
                subset = popDF[popDF.fitness > self.min_fitness]
                # Drop columns with all NA values before concatenation
                self.data = pd.concat([self.data, subset], ignore_index=True)
                self._update_data_frame()
                if self.data.shape[0] > self.max_size:
                    self.data = self.data[:self.max_size]

    def _update_data_frame(self):
        """
        Sort, drop duplicates, and update the minimum fitness in the hall of
        fame data.
        """
        self.data = self.data.sort_values("fitness", ascending=False)
        self.data = self.data.drop_duplicates(keep="first", ignore_index=True)
        self.data = self.data.reset_index(drop=True)
        self.custom_drop()
        self.min_fitness = self.data["fitness"].min()

    def custom_drop(self):
        """
        Perform custom operations on the hall of fame data.

        This method applies custom transformations and drops duplicates in the
        data. For each variable, it updates the weight, trans, and shape
        columns by multiplying them with the variable's value and performs
        additional adjustments.
        """
        for v in self.variables:
            v_str = str(v)
            self.data[f"{v_str}_weight"] = (
                self.data[v_str] * self.data[f"{v_str}_weight"]
            )
            self.data[f"{v_str}_trans"] = (
                self.data[v_str] * self.data[f"{v_str}_trans"]
            )
            self.data[f"{v_str}_shape"] = (
                self.data[v_str] * self.data[f"{v_str}_shape"]
            )
            self.data[f"{v_str}_asym"] = (
                self.data[v_str] * self.data[f"{v_str}_asym"]
            )
            temp = self.data[f"{v_str}_trans"]
            temp[temp > 1] = 1
            self.data[f"{v_str}_shape"] = self.data[v_str] * temp

        self.data = self.data.drop_duplicates(keep="first", ignore_index=True)
        self.data = self.data.reset_index(drop=True)

    def printHOF(self, max_row=None, max_col=None):
        """
        Print the hall of fame data sorted by fitness.

        Args:
            max_row (int, optional): Maximum number of rows to display.
                                     Defaults to None.
            max_col (int, optional): Maximum number of columns to display.
                                     Defaults to None.
        """
        self.data = self.data.sort_values("fitness", ascending=False)
        self.data = self.data.reset_index(drop=True)
        with pd.option_context("display.max_rows",
                               max_row,
                               "display.max_columns",
                               max_col):
            print(self.data)

    def printBest(self, max_row=None, max_col=None):
        """
        Print the best model in the hall of fame.

        Args:
            max_row (int, optional): Maximum number of rows to display.
                                     Defaults to None.
            max_col (int, optional): Maximum number of columns to display.
                                     Defaults to None.
        """
        if self.best is None:
            self.get_best_model()
        self.best = self.best.sort_values("Variable", ascending=False)
        self.best = self.best.reset_index(drop=True)
        with pd.option_context("display.max_rows",
                               max_row,
                               "display.max_columns",
                               max_col):
            print(self.best)

    def printRVI(self, max_row=None, max_col=None):
        """
        Print the relative variable importance (RVI) for models in the hall of
        fame.

        Args:
            max_row (int, optional): Maximum number of rows to display.
                                     Defaults to None.
            max_col (int, optional): Maximum number of columns to display.
                                     Defaults to None.
        """
        if self.rvi is None:
            self.relative_variable_importance()
        self.rvi = self.rvi.sort_values("RVI", ascending=False)
        self.rvi = self.rvi.reset_index(drop=True)
        with pd.option_context("display.max_rows",
                               max_row,
                               "display.max_columns",
                               max_col):
            print(self.rvi)

    def printMAW(self, max_row=None, max_col=None):
        """
        Print the model average weights (MAW) for models in the hall of fame.

        Args:
            max_row (int, optional): Maximum number of rows to display.
                                     Defaults to None.
            max_col (int, optional): Maximum number of columns to display.
                                     Defaults to None.
        """
        if self.maw is None:
            self.model_average_weights()
        self.maw = self.maw.sort_values("MAW", ascending=False)
        self.maw = self.maw.reset_index(drop=True)
        with pd.option_context("display.max_rows",
                               max_row,
                               "display.max_columns",
                               max_col):
            print(self.maw)

    def recalculate_aic(self):
        """
        Recalculate the Akaike Information Criterion (AIC) for the models.

        AIC is calculated using the formula: -2 * log-likelihood + 2 * k,
        where k is the number of variables in the model.
        """
        k = self.data[self.variables].apply(lambda x: (x == 1).sum(), axis=1)
        self.data["aic"] = -2 * self.data["loglik"] + 2 * k

    def calculate_bic(self, n):
        """
        Calculate the Bayesian Information Criterion (BIC) for the models.

        Args:
            n (int): The number of observations.

        BIC is calculated using the formula: -2 * log-likelihood + k * log(n),
        where k is the number of variables in the model.
        """
        k = self.data[self.variables].apply(lambda x: (x == 1).sum(), axis=1)
        self.data["bic"] = -2 * self.data["loglik"] + k * np.log(n)

    def delta_aic(self):
        """
        Calculate the difference in AIC values relative to the best model.

        This method calculates the delta AIC, which is the difference between
        each model's AIC and the lowest AIC in the dataset.
        """
        if self.data.shape[0] <= 0:
            return

        if "delta_aic_best" in self.data.columns:
            return

        self.data["aic"] = self.data["aic"] * -1  # Reverse the neg sign
        best_aic = self.data["aic"].min()
        self.data["delta_aic_best"] = self.data["aic"] - best_aic

    def akaike_weights(self):
        """
        Calculate and assign Akaike weights to the models in the hall of fame.

        This method computes the Akaike weights based on the delta AIC values
        of the models.

        The weight of each model is calculated using the formula:
        e^(-1/2 * delta_aic_best) / sum(e^(-1/2 * delta_aic_best(k)))
        where the denominator is summed over all models.
        """
        if self.data.shape[0] <= 0:
            return

        if "delta_aic_best" not in self.data.columns:
            self.delta_aic()

        exp_delta_aic = np.exp(-0.5 * self.data["delta_aic_best"])
        self.data["akaike_weight"] = exp_delta_aic / exp_delta_aic.sum()

    def cumulative_akaike(self, threshold=1.0):
        """
        Calculate the cumulative Akaike weights and identify models to keep
        based on a threshold.

        This method calculates the cumulative sum of Akaike weights and marks
        models to keep based on the specified cumulative weight threshold.

        Args:
            threshold (float, optional): The threshold for cumulative Akaike
                                         weight. Defaults to 1.0.
        """
        if self.data.shape[0] <= 0:
            return

        self.data = self.data.reset_index(drop=True)
        threshold = float(threshold)

        if "akaike_weight" not in self.data.columns:
            self.akaike_weights()

        self.data["acc_akaike_weight"] = self.data["akaike_weight"].cumsum()

        if 0.0 < threshold < 1.0:
            max_weight = self.data["acc_akaike_weight"].max()
            if max_weight > threshold:
                cutoff = self.data[
                        self.data["acc_akaike_weight"].gt(threshold)
                    ].index[0]
                self.data["keep"] = ["True"] * (cutoff + 1) + \
                    ["False"] * (self.data.shape[0] - (cutoff + 1))
            else:
                self.data["keep"] = ["True"] * self.data.shape[0]
        else:
            self.data["keep"] = ["True"] * self.data.shape[0]

    def relative_variable_importance(self, ignore_keep=False):
        """
        Calculate the relative variable importance (RVI) based on Akaike
        weights.

        This method computes the sum of weighted Akaike weights for each
        variable and ranks them, providing an insight into the relative
        importance of each variable.

        Args:
            ignore_keep (bool, optional): If True, includes all models in the
                                          calculation regardless of their
                                          'keep' status. Defaults to False.
        """
        # Clear previous calculations
        self.rvi = pd.DataFrame(columns=["variable", "RVI"])

        # Selecting the subset of data based on 'keep' column
        sub = self.data if ignore_keep else self.data[self.data.keep == "True"]

        # Prepare data for concatenation
        rows = []
        for v in self.variables:
            sum_weights = (sub[v] * sub["akaike_weight"]).sum()
            rows.append({"variable": v, "RVI": sum_weights})

        # Concatenate all rows into the DataFrame
        self.rvi = pd.concat(
            [self.rvi, pd.DataFrame(rows)], ignore_index=True
        )

        # Sorting and resetting index
        self.rvi = self.rvi.sort_values(
            "RVI", ascending=False
        ).reset_index(drop=True)

    def get_best_model(self):
        """
        Retrieve the best model based on fitness.

        This method identifies the best model in the hall of fame and extracts
        its variable weights, shapes, and transformations.

        Returns:
            pd.DataFrame: A DataFrame representing the best model with its
                          variables and attributes.
        """
        if self.data.empty:
            return None

        best_model = self.data.iloc[0]
        best_model_variables = [
            {
                "Variable": var,
                "Weight": best_model[f"{var}_weight"],
                "Shape": best_model[f"{var}_shape"],
                "Transformation": best_model[f"{var}_trans"],
                "Asymmetry": best_model[f"{var}_asym"],
            }
            for var in self.variables
        ]

        self.best = pd.DataFrame(best_model_variables)
        return self.best

    def model_average_weights(self, ignore_keep=False):
        """
        Calculate the model average weights (MAW) based on Akaike weights.

        This method computes the sum of weighted Akaike weights for each
        variable and ranks them, providing an insight into the average weight
        of each variable across models.

        Args:
            ignore_keep (bool, optional): If True, includes all models in the
                                          calculation regardless of their
                                          'keep' status. Defaults to False.
        """
        # Clear previous calculations
        self.maw = pd.DataFrame(columns=["variable", "MAW"])

        # Selecting the subset of data based on 'keep' column or ignoring it
        sub = self.data if ignore_keep else self.data[self.data.keep == "True"]

        # Prepare data for concatenation
        rows = []
        for v in self.variables:
            weight_col = f"{v}_weight"
            sum_weights = (sub[weight_col] * sub["akaike_weight"]).sum()
            rows.append({"variable": v, "MAW": sum_weights})

        # Concatenate all rows into the DataFrame
        self.maw = pd.concat([self.maw, pd.DataFrame(rows)], ignore_index=True)

        # Sorting and resetting index
        self.maw = self.maw.sort_values(
            "MAW", ascending=False
        ).reset_index(drop=True)

    def cleanHOF(self):
        """
        Clean and prepare the hall of fame data for analysis or presentation.

        This method sorts the data by fitness and resets the index. For each
        variable, it replaces weight, transformation, and shape values with NaN
        where the variable value is 0.

        Returns:
            pd.DataFrame: The cleaned hall of fame data.
        """
        ret = self.data.sort_values("fitness", ascending=False)
        ret = ret.reset_index(drop=True)

        for v in self.variables:
            mask = ret[v] == 0
            for attr in ["_weight", "_trans", "_shape", "_asym"]:
                ret[f"{v}{attr}"][mask] = np.nan

        if ret.iloc[0]["fitness"] == (ret.iloc[0]["aic"] * -1):
            ret["fitness"] = ret["fitness"] * -1

        return ret

    def get_variables(self):
        """
        Get the list of variables used in the hall of fame models.

        Returns:
            list: A list of variables.
        """
        return self.variables

    def getBest(self):
        """
        Retrieve the best model data.

        If not already calculated, this method first calculates the best model.
        It then sorts and resets the index of the best model data.

        Returns:
            pd.DataFrame: DataFrame containing the best model data.
        """
        if self.best is None:
            self.get_best_model()
        self.best = self.best.sort_values(
            "Variable", ascending=False).reset_index(drop=True)
        return self.best

    def getRVI(self):
        """
        Retrieve the relative variable importance (RVI) data.

        If not already calculated, this method first calculates the RVI.
        It then sorts and resets the index of the RVI data.

        Returns:
            pd.DataFrame: DataFrame containing the RVI data.
        """
        if self.rvi is None:
            self.relative_variable_importance()
        self.rvi = self.rvi.sort_values(
            "RVI", ascending=False).reset_index(drop=True)
        return self.rvi

    def getMAW(self):
        """
        Retrieve the model average weights (MAW) data.

        If not already calculated, this method first calculates the MAW.
        It then sorts and resets the index of the MAW data.

        Returns:
            pd.DataFrame: DataFrame containing the MAW data.
        """
        if self.maw is None:
            self.model_average_weights()
        self.maw = self.maw.sort_values(
            "MAW", ascending=False).reset_index(drop=True)
        return self.maw

    def getHOF(self, only_keep=False):
        """
        Retrieve the Hall of Fame data.

        This method sorts the Hall of Fame data by fitness and resets the index
        If the only_keep flag is set to True, it filters the data to return
        only the models marked as 'keep'.

        Args:
            only_keep (bool, optional): Flag to return only the models marked
                                        as 'keep'. Defaults to False.

        Returns:
            pd.DataFrame: DataFrame containing the Hall of Fame data, filtered
                          if only_keep is True.
        """
        self.data = self.data.sort_values(
            "fitness", ascending=False).reset_index(drop=True)
        return self.data[self.data.keep == "True"] if only_keep else self.data

    def plot_ICprofile(self, oname="out", diff=2):
        """
        Create and save a scatter plot of the information criterion profile.

        The plot displays models ordered by AIC values, with a horizontal line
        representing the specified AIC difference threshold.

        Args:
            oname (str, optional): The output filename prefix for the saved
                                   plot. Defaults to "out".
            diff (int, optional): The AIC difference threshold for highlighting
                                  the plot. Defaults to 2.
        """
        diff = int(diff)
        # Sort data by AIC and create model numbering
        dat = self.data.sort_values("aic", ascending=True).round(3)
        dat.reset_index(inplace=True, drop=True)
        dat["model"] = dat.index + 1

        # Plotting
        sns.set(style="ticks")
        p = sns.scatterplot(
            data=dat, x="model", y="aic", hue="r2m", size="r2m",
            style="keep", alpha=0.6
        )
        p.axhline(dat["aic"].min() + diff, ls="--", c="red")
        plt.title("IC Profile")

        # Save and clear the plot
        plt.savefig(f"{oname}.ICprofile.pdf")
        plt.clf()
        plt.close()

    def plotMetricPW(self, oname="out"):
        """
        Create and save a pair plot of various metrics.

        The plot visualizes pairwise relationships in the dataset for selected
        columns. If the 'akaike_weight' column exists, it is included in the
        plot.

        Args:
            oname (str, optional): The output filename prefix for the saved
                                   plot. Defaults to "out".
        """
        # Determine columns to plot
        cols = ["aic", "loglik", "r2m", "delta_aic_null", "keep"]
        if "akaike_weight" in self.data.columns:
            cols.append("akaike_weight")

        # Subset the data for plotting
        dat = self.data[cols]

        # Plotting
        sns.set(style="ticks")
        sns.pairplot(dat, hue="keep", kind="scatter")
        plt.title("Pair plot")

        # Save and clear the plot
        plt.savefig(f"{oname}.pairPlot.pdf")
        plt.clf()
        plt.close()

    def plotVariableImportance(self, oname="out", cutoff=0.8):
        """
        Create and save a bar plot of variable importance.

        The plot displays the relative variable importance (RVI) as a bar plot
        with a cutoff line indicating the specified threshold.

        Args:
            oname (str, optional): The output filename prefix for the saved
                                   plot. Defaults to "out".
            cutoff (float, optional): The threshold for highlighting important
                                      variables. Defaults to 0.8.
        """
        cutoff = float(cutoff)

        # Prepare data for plotting
        sub = self.rvi.sort_values("RVI", ascending=False)

        # Plotting
        sns.set(style="ticks")
        p = sns.barplot(data=sub, x="RVI", y="variable")
        p.axvline(cutoff, ls="--", c="red")
        plt.title("Importance of Terms")

        # Save and clear the plot
        plt.savefig(f"{oname}.varImportance.pdf")
        plt.clf()
        plt.close()

    def plotModelAveragedWeights(self, oname="out"):
        """
        Create and save a bar plot of model-averaged weights.

        The plot displays the model-averaged weights (MAW) as a bar plot,
        illustrating the average weight of each variable across all models.

        Args:
            oname (str, optional): The output filename prefix for the saved
                                   plot. Defaults to "out".
        """
        # Prepare data for plotting
        sub = self.maw.sort_values("MAW", ascending=False)

        # Plotting
        sns.set(style="ticks")
        sns.barplot(data=sub, x="MAW", y="variable")
        plt.title("Model-Averaged Weights")

        # Save and clear the plot
        plt.savefig(f"{oname}.modavgWeights.pdf")
        plt.clf()
        plt.close()

    def writeModelSummary(self, oname):
        """
        Write the Hall of Fame model summary to a TSV file.

        Args:
            oname (str): The output filename prefix for the saved TSV file.
        """
        out_df = self.cleanHOF()
        out_df.to_csv(
            f"{oname}.HallOfFame.tsv", sep="\t", index=False, na_rep="-"
        )

    def writeMAW(self, oname):
        """
        Write the model-averaged weights (MAW) to a TSV file.

        Args:
            oname (str): The output filename prefix for the saved TSV file.
        """
        if self.maw is None:
            self.model_average_weights()
        self.maw.to_csv(
            f"{oname}.modavgWeights.tsv", sep="\t", index=False, na_rep="-"
        )

    def writeRVI(self, oname):
        """
        Write the relative variable importance (RVI) to a TSV file.

        Args:
            oname (str): The output filename prefix for the saved TSV file.
        """
        if self.rvi is None:
            self.relative_variable_importance()
        self.rvi.to_csv(
            f"{oname}.varImportance.tsv", sep="\t", index=False, na_rep="-"
        )

    def writeBest(self, oname):
        """
        Write the best model details to a TSV file.

        Args:
            oname (str): The output filename prefix for the saved TSV file.
        """
        if self.best is None:
            self.get_best_model()
        self.best.to_csv(
            f"{oname}.bestModel.tsv", sep="\t", index=False, na_rep="-"
        )

    @classmethod
    def from_dataframe(cls, df, max_size=None):
        """
        Create a HallOfFame instance from a DataFrame.

        Args:
            cls (class): The class constructor.
            df (pd.DataFrame): The DataFrame containing the data.
            max_size (int, optional): The maximum size for the HallOfFame.
                                      Defaults to None.

        Returns:
            HallOfFame: A new HallOfFame instance initialized with data from
                        the DataFrame.
        """
        # Extract variable names from the DataFrame column names
        variable_prefixes = ["_weight", "_trans", "_shape", "_asym"]
        variables = set()
        for col in df.columns:
            for prefix in variable_prefixes:
                if col.endswith(prefix):
                    variables.add(col.replace(prefix, ""))
        variables = list(variables)
        if "acc_akaike" in variables:
            variables.remove("acc_akaike")
        if "akaike" in variables:
            variables.remove("akaike")

        df = df.fillna(0.0)
        # Set max_size to the maximum allowable row number if not provided
        if max_size is None:
            max_size = pd.options.display.max_rows

        # Create a new instance of the HallOfFame class
        instance = cls(variables, max_size)

        df["aic"] = df["aic"] * -1  # Reverse the neg sign added for maximizing
        best = df["aic"].min()
        df["delta_aic_best"] = df["aic"] - best

        # Update the instance attributes with the provided DataFrame
        instance.data = df.copy()

        return instance

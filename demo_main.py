import matplotlib.pyplot as plt
from geo_utils.combine import combine_and_crop


if __name__ == "__main__":
    s57_candidates_2 = combine_and_crop(
        lati=34.754640848405174, longi=119.37176328480543, size=2000,
    )
    s57_candidates_2.plot()
    plt.show()
I ran the script with a direction of negative I believe and an N of 1000 for the null distribution (which takes a while) and a Threshold of 100.

A threshold of 100 is arguably not the right thing to do, but it did not seem to change the scores/stats that much whether running with 10 or 100.

The correct way to do it now, I think, is to run the analysis with a threshold of 50~60 for both directions, so you do the analysis twice.
Once for a ranking of positive to negative and once for a ranking of negative to positive, both with a threshold of 50-60. 
Combining both ranking analysis will then be the final result, in which you take the q-value for the right direction for a given gene (e.g. a depleted gene has the negative direction q-value while a postive has the positive direction q-value).
# Conflict-Task-Analysis-Package
## Task Description
Photometry analysis for the platform mediated conflict task, in which mice must balance their desire for a reward with the risk of electric shock (behavioral task adapted for mice by Christian Bravo). In this task, mice learn to avoid a mild shock by taking shelter on top of an insulated platform. A 20 second beep plays directly before the shock, giving mice time to retreat to the platform. In the low-conflict version of this task, a water or sugar water reward is available at a water port opposite the platform at all times, so the mice must leave the safety of the platform for the reward, but can do so without much risk (since they know the shock won't start unless the sound plays). In the high conflict task, the liquid reward is only available during the 20 second tone predicting the shock. This requires the mice to balance the increased risk with their desire for the reward. 

<p align="center">
  <img width=800 src="https://github.com/SaraEBoyle/Conflict-Task-Analysis-Package/assets/83416542/964ee017-cfc1-4dba-bc8e-f8b7d6d0917c" alt="Task Diagram" />
</p>

## Analysis Summary

My analysis package takes video from recorded sessions of this task, tracks the location of the mouse, and calculates the time mice spend in the safety of the platform vs. in front of the water reward. By tracking the location of the mice during different timepoints, we can track how individual mice respond to risky situations. Here's an example of the tracking in a high conflict session in which the mouse receives a shock:

<p align="center">
  <img width=800 src="https://github.com/SaraEBoyle/Conflict-Task-Analysis-Package/assets/83416542/7bd1aa71-575d-4b5c-87bd-b38bcd418b22" alt="mouse_trackingmouse_tracking" />
</p>

## Tracking the Mice in Low and High Conflict

For each individual mouse, my code tracks the location through the entire session and plots the trajectory of the mouse. I also plot heatmaps showing the areas where the mice spend the most time. My code splits the tracking data into periods with no sound playing (no shock risk) and periods with sound playing (imminent shock risk). In this example mouse, you can see in the low conflict (LC) session that the mouse spends the majority of its time on the platform when the sound is playing (which makes sense, since the platform is safe). However, when the sound stops, the mouse leaves the platform to drink from the reward spout. This can also be plotted as the probability that the mouse is in the platform zone at any given time point (calculated using all recorded trials), and predictably, p(Platform) increases at the onset of the sound. 

In high conflict (HC) sessions, the reward is only available during the sound period. So, this particular mouse chooses to spend more time in front of the reward port during the sound, a much riskier behavior than in the low conflict session. This allows the mouse to receive rewards, but resulted in it receiving a shock in one of the trials. This behavior is also reflected in the number of times this mouse poked its nose into the reward port (which I measured using a beam breaker). 

![platform likelihood and tracking](https://github.com/SaraEBoyle/Conflict-Task-Analysis-Package/assets/83416542/1b9cb4ec-0b4d-473c-89e4-c69e5deccb0c)


## Using the Conflict Task to Understand Risky Behavior

By training many mice on this task, we can start to understand how different mice juggle risk and reward, and we hope to understand the underlying brain circuits that make some mice risk-prone, and others risk-averse. We hope to use that knowledge to learn more about how humans might engage in risky decisions, like gambling or drug use.

Here is an average plot of all of the mice I trained on this task. All mice showed some of the same basic behaviors as the example mouse: during the low conflict task, they tend to go to the platform as soon as the sound starts, since there is no incentive for them to risk a shock. However, during high conflict sessions, the mice tend to stay away from the platform during the sound longer in order to get a reward. 

![group behavior alone](https://github.com/SaraEBoyle/Conflict-Task-Analysis-Package/assets/83416542/f13f205b-997c-467d-83e9-b0a8c691fdba)

## Individual Variation in Willingness to Take Risks

Although most mice behave the same during low conflict sessions, there is a lot of variation across individuals during the high conflict task. We can see this by plotting the percent time spent in the reward zone vs. the percent time spent in the platform zone. During low conlfict, when the sound is off (no risk) mice spend most of their time in the reward zone, and little time in the platform zone. When the sound begins to play, this reverses, as all mice prefer the safety of the platform.

During high conflict sessions, when the sound is off, the mice can't get the reward, though many of them choose to wait in the reward zone. However, once the sound begins to play, we see a much larger variation in strategies. Some mice are much more risk prone (RP), choosing to spend as much time as possible in the reward zone, even if it means getting shocked more often. Some mice are risk-averse (RA), almost immediately returning to the platform once the sound plays, and receiving very little reward as a result.

![risk averse vs risk prone alone](https://github.com/SaraEBoyle/Conflict-Task-Analysis-Package/assets/83416542/cc98c835-bd60-4990-8b77-f8b478cc9000)

## Differences in Brain Activity Between Risk-Prone and Risk-Averse Mice

We are beginning to see some differences in the brains of risk-averse vs. risk-prone mice. This analysis package can also handle photometry recordings (a measure of bulk brain activity in a specific brain area). We've found an area that responds more strongly to the sound in risk-averse mice than in risk-prone mice (can't tell you which brain region, it's a secret, shhh).

With this analysis package, we can compare brain responses to major events in the conflict task, like the sound onset (CS, in these graphs) and the shock. We can see that this brain area doesn't have any difference in shock response across the high conflict and the low conflict sessions. However, the risk-averse mice show sound responses only in the high conflict sessions. This could indicate that this brain area biases mice to avoid risks, or weigh avoiding a risk as more important than seeking a reward. 

![photometry responses alone](https://github.com/SaraEBoyle/Conflict-Task-Analysis-Package/assets/83416542/d383d711-a989-434c-86d9-8099a9e0031f)

## How to Use This Package

### Inputs:
1. Analog input file (.csv)- file containing timing data for syncing each bpod trial to photometry data. Ensures timing is stable throughout entire session.
2. Photometry data file (.csv)- raw data array (n x m) where n is number of data entries, m = 1 is time, and m = 2, 3, ... is fluorescence data from each fiber
3. Bpod file (.mat)- contains all event and state data from bpod events
4. Centroid file (.csv)- Tracking data obtained in real time using Bonsai workflow (included)
5. Background picture of arena (.jpeg) before placing the mouse inside, for use in background subtraction for better movement tracking

### Outputs:
1. Probability of each mouse to enter the platform region, across all trials
2. Risk index (% time in reward zone vs. % time on platform) of the session, split between sound on and sound off time periods
3. Tracking data traces and heatmaps
4. Aligned photometry responses to all major events (sound onset, shock, platform entry, platform exit, reward received, etc.)

### Instructions:
1. Complete conflict task session (photometry recording optional). Use included Bpod behavioral protocol for experiment control and Bonsai file for data acquisition. See Bpod file for experiment assembly instructions/camera and Bonsai setup.
2. After acquiring the data, run full_photometry_conflict_analysis.m from folder containing required inputs. There are a variety of parameters described that you can adjust.
3. If recording photometry responses, follow the prompts for optionally trimming the session/completeing bleaching correction.





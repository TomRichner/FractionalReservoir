% Independent synthetic checks: noise, localized events, and global burst.
t=0:.005:10; x=.01*randn(100,numel(t)); r=zeros(size(x));
groups=repelem((1:5)',20);
quiet=ied_detect_events(t,x,r,groups,struct());
assert(quiet.metrics.event_count==0);
x(1:20,:)=x(1:20,:)+.6*exp(-((t-3)/.09).^2);
x(41:60,:)=x(41:60,:)+.6*exp(-((t-6)/.09).^2);
localized=ied_detect_events(t,x,r,groups,struct());
assert(localized.metrics.event_count==2);
assert(all(abs([localized.events.time]-[3 6])<.10));
assert(localized.metrics.global_fraction==0);
x=x+.6*exp(-((t-8)/.09).^2);
global_run=ied_detect_events(t,x,r,groups,struct());
assert(global_run.metrics.event_count==3 && global_run.metrics.global_fraction>0);
disp('IED candidate detector: noise, localized events, global burst checks passed.');

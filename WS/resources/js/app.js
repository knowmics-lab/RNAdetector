require('./bootstrap');

window.$(function () {
    window.$('.logout-button').click(function (e) {
        e.preventDefault();
        window.$('#logout-form').submit();
    });
});
